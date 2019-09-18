#include "vars.h"
#include "proto.h"

int main(int argc, char *argv[]) {

  // Set up MPI
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  fftw_mpi_init();

  const int bsz = 2000;
  char buf[bsz], fout_name[bsz], fout_name_2D[bsz], fin_name[bsz];

  if(argc < 5) {
    if (ThisTask == 0) {
      printf("\nParameters are missing.\n");
      printf("Call with <InFile, OutFile, do_2D, do_log>\n\n");
    }
    FatalError("read_data", 55);
  }

  sprintf(fin_name,  "%s", argv[1]);
  sprintf(fout_name, "%s", argv[2]);
  do_2D = atoi(argv[3]);
  do_log = atoi(argv[4]);

  if (ThisTask == 0) {
    printf("Running power spectrum with do_2D=%d, do_log=%d\n", do_2D, do_log); fflush(stdout);
    if (do_log) {
      if (mink == 0.0) { printf("Error mink=0 with log-binning\n");  FatalError("read_data", 64); }
    }
    if ((do_interlacing != 0) && (do_interlacing != 1)) { printf("ERROR: do_interlacing must be 0 or 1\n");  FatalError("read_data", 71); }
    if ((interp_order < 1) && (interp_order > 4)) { printf("ERROR: interp_order must be between 1 and 4\n");  FatalError("read_data", 72); }
    printf("Allocating grid...\n"); 
    fflush(stdout);
  }

  // *********************************************************  
  const double dx  = (XMAX-XMIN)/(double)NX;
  const double dy  = (YMAX-YMIN)/(double)NY;
  const double dz  = (ZMAX-ZMIN)/(double)NZ;
  ptrdiff_t ii, jj, kk;

  ptrdiff_t Local_nx;       // The number of slices on the task
  ptrdiff_t Local_x_start;  // The global start of the slices on the task
  ptrdiff_t alloc_local = fftw_mpi_local_size_3d(NX, NY, NZ/2+1, MPI_COMM_WORLD, &Local_nx, &Local_x_start);
  ptrdiff_t NTOT = 2*alloc_local;
  // Add additional planes as needed based on the interpolation order and whether or not we are interlacing
  ptrdiff_t alloc_slice = 2*NY*(NZ/2+1);
  ptrdiff_t last_slice = Local_nx*alloc_slice;
  //ptrdiff_t nextra = interp_order;
  ptrdiff_t Local_nxtra = Local_nx+interp_order;
  Total_size = NTOT+interp_order*alloc_slice;

  // Set the neighbouring tasks
  int * Local_nx_table = (int *)malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);
  if (Local_nx == 0) {
    LeftTask = MPI_PROC_NULL;
    RightTask = MPI_PROC_NULL;
  } else {
    LeftTask = ThisTask;
    do {
      LeftTask--;
      if(LeftTask < 0) LeftTask = NTask - 1;
    } while(Local_nx_table[LeftTask] == 0);
      
    RightTask = ThisTask;
    do {
      RightTask++;
      if(RightTask >= NTask) RightTask = 0;
    } while(Local_nx_table[RightTask] == 0);
  }

  double * ddg = (double*)calloc(Total_size,sizeof(double));
  fftw_plan dp_r2c = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg,(fftw_complex*)ddg,MPI_COMM_WORLD,FFTW_ESTIMATE); 

  double * ddg_2;
  fftw_plan dp_r2c_2;
  if (do_interlacing) {
    ddg_2 = (double*)calloc(Total_size,sizeof(double));
    dp_r2c_2 = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_2,(fftw_complex*)ddg_2,MPI_COMM_WORLD,FFTW_ESTIMATE);   
  }

  // *********************************************************
  // Normally this is a quick parallel read, but this is unnecessary for the MockChallenge
  // (there is only one small file), instead, every processor just reads in the file 
  // simultaneously and only saves galaxies that are in the correct portion of the grid.

  if (ThisTask == 0) {printf("Reading data...\n"); fflush(stdout);}

  FILE * fp;
  double NIN=0, NPART_TOT=0;
  double XMIN_LOCAL = Local_x_start*dx+XMIN;
  double XMAX_LOCAL = (Local_x_start+Local_nx)*dx+XMIN;
  if((fp=fopen(fin_name,"r"))==NULL) { printf("Task %d cannot open input file\n", ThisTask);  FatalError("read_data", 96); }
  while(fgets(buf,bsz,fp)) {

    if(strncmp(buf,"#",1)==0) continue;

    double tx,ty,tz,tz_rsd;
    if(sscanf(buf,"%lf %lf %lf %lf\n",&tx,&ty,&tz,&tz_rsd)!=4) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
    tz = tz_rsd;

    double w = 1.0;

    NPART_TOT += w;

    if ((tx < XMIN) || (tx >= XMAX) || (ty < YMIN) || (ty >= YMAX) || (tz < ZMIN) || (tz >= ZMAX)) {
      printf("out of bounds: %lf %lf %lf\n", tx, ty, tz);
      FatalError("read_data", 110);
    }

    NIN += add_to_grid(tx, ty, tz, w, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, dx, dy, dz, ddg);
    if (do_interlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, dx, dy, dz, ddg_2);

  }
  printf("Task %d read %lf particles, gridded %lf\n", ThisTask, NPART_TOT, NIN);
  fflush(stdout);

  // Copy across the extra slices from the task on the left and add it to the leftmost slices
  // of the task on the right. Skip over tasks without any slices.
  if (interp_order > 0) {
    MPI_Status status;
    double * temp_ddg = (double *)calloc(interp_order*alloc_slice,sizeof(double));
    ierr = MPI_Sendrecv(&(ddg[last_slice]),interp_order*alloc_slice*sizeof(double),MPI_BYTE,RightTask,0,
                        &(temp_ddg[0]),interp_order*alloc_slice*sizeof(double),MPI_BYTE,LeftTask,0,MPI_COMM_WORLD,&status);
    for (int i=0;i<interp_order*alloc_slice;i++) ddg[i] += temp_ddg[i];
    if (do_interlacing) {
      ierr = MPI_Sendrecv(&(ddg_2[last_slice]),interp_order*alloc_slice*sizeof(double),MPI_BYTE,RightTask,0,
                          &(temp_ddg[0]),interp_order*alloc_slice*sizeof(double),MPI_BYTE,LeftTask,0,MPI_COMM_WORLD,&status);
      for (int i=0;i<interp_order*alloc_slice;i++) ddg_2[i] += temp_ddg[i];
    }
    free(temp_ddg);
  }

  // *********************************************************
  FILE *fout, *fout_2D;
  if (ThisTask == 0) {
    if((fout=fopen(fout_name,"w"))==NULL) { printf("cannot open output file: %s\n", fout_name); FatalError("read_data", 142); }
    printf("%s\n",fout_name);
    if (do_2D) {
      sprintf(fout_name_2D, "%s_2D", fout_name);
      if((fout_2D=fopen(fout_name_2D,"w"))==NULL) { printf("cannot open output file: %s\n", fout_name_2D); FatalError("read_data", 142); }
      printf("%s\n",fout_name_2D);
    }
  }

  // subtract off mean density
  // Referenced to global mean
  double NIN_glob = 0;
  MPI_Allreduce(&NIN, &NIN_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double density_grid = NIN_glob/((double)NX*(double)NY*(double)NZ);
  double nbar = NIN_glob/((XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN));
  double nbwsq = NIN_glob;
  double nbsqwsq = nbar*NIN_glob;
  
  double sumddg = 0.0;
  for (ii=0; ii<Local_nx; ii++) {
    for (jj=0; jj<NY; jj++) {
      for (kk=0; kk<NZ; kk++) {
        long long ind = kk+2*(NZ/2+1)*(jj+NY*ii);
        sumddg += ddg[ind];
        ddg[ind] -= density_grid;
        if (do_interlacing) ddg_2[ind] -= density_grid;
      }
    }
  }
  printf("%lf\n", sumddg);

  if (ThisTask == 0) {
    printf("# Nparticles = %g\n", NIN_glob);
    printf("# integrals: nbar    = %g\n",nbar);
    printf("# integrals: nbwsq   = %g\n",nbwsq);
    printf("# integrals: nbsqwsq = %g\n",nbsqwsq);
    fflush(stdout);
    fprintf(fout, "# k, pk0, pk2, pk4, Nmodes\n");
    if (do_2D) {
      fprintf(fout_2D, "# k, mu, pk, Nmodes\n");
    }
  }
     
  // *********************************************************
  // Set up bins for P(k) data 
  double binwidth;
  if (do_log) {
    binwidth=(log10(maxk)-log10(mink))/(double)NB; // bin width
  } else {
    binwidth=(maxk-mink)/(double)NB; // bin width
  }

  double * Pg0_2D, * Nk_2D;
  if (do_2D) { 
    Pg0_2D = (double*)calloc(NB*NMU, sizeof(double));   // binned power spectrum
    Nk_2D = (double*)calloc(NB*NMU, sizeof(double));        // number of modes in each bin
  }
  double * Pg0 = (double*)calloc(NB, sizeof(double));       // binned power spectrum
  double * Pg2 = (double*)calloc(NB, sizeof(double));       // binned power spectrum
  double * Pg4 = (double*)calloc(NB, sizeof(double));       // binned power spectrum
  double * Nk = (double*)calloc(NB, sizeof(double));        // number of modes in each bin
  int * Pgfill = (int*)calloc(NB, sizeof(int));          // check if P(k) bin already filled
  
  // calculate Nyquist frequency
  double ny_x = (M_PI/dx);
  double ny_y = (M_PI/dy);
  double ny_z = (M_PI/dz);
  double min_nyquist=ny_x;
  if(ny_y<min_nyquist) min_nyquist=ny_y;
  if(ny_z<min_nyquist) min_nyquist=ny_z;

  // Specify the LOS
  double vec[3] = {0.0, 0.0, 1.0};

  if (ThisTask == 0) {
    printf("Now Fourier transforming overdensity field\n");
    fflush(stdout);
  }
  fftw_execute(dp_r2c);
  if (do_interlacing) fftw_execute(dp_r2c_2);

  // Work out |k| for each frequency component
  if (ThisTask == 0) {
    printf("Iterating over grid cells\n");
    fflush(stdout);
  }

  double sx = 1.0/((float)NX*dx), sy = 1.0/((float)NY*dy), sz = 1.0/((float)NZ*dz)
  for (int i=Local_x_start; i<Local_x_start+Local_nx; i++) {
    for (int j=0; j<NY; j++) {
      for (int k=0; k<=NZ/2; k++) {
	
        // frequency in x
        double fx;
        if (i <= NX/2) { 
          fx = (float)i*sx;
        } else {
          fx = ((float)i-(float)NX)*sx;
        }

        // frequency in y
        double fy;
        if (j <= NY/2) {
          fy = (float)j*sy;
        } else {
          fy = ((float)j-(float)NY)*sy;
        }

	      // frequency in z
	      double fz = (float)k*sz;
	
        // length of k vector, value of mu and bin
        double fktot=2.0*M_PI*sqrt(fx*fx+fy*fy+fz*fz);
        double mu;
        if (fktot > 0) {
          if (los == 1) {
            mu = 2.0*M_PI*fx/fktot;
          } else if (los == 2) {
            mu = 2.0*M_PI*fy/fktot;
          } else {
            mu = 2.0*M_PI*fz/fktot;
          }
        } else {
          mu = 0.0;
        }
        long long kbin;
        if (do_log) {
          kbin = (long long)((float)(log10(fktot)-log10(mink))/binwidth);
        } else {
          kbin = (long long)((float)(fktot-mink)/binwidth);
        }
        long long mubin = (long long)(mu*(float)NMU);
        long long bin = kbin*NMU + mubin;

	      if ((kbin >= 0) && (kbin < NB) && (fktot < min_nyquist) && (Pgfill[kbin]==0)) {
	  
	  	    // set up correction for gridding - in effect we're
	        // convolving the density field with a top-hat function in
	        // each direction, so we're multiplying each Fourier mode by
	        // a sinc function. To correct this, we therefore divide by
	        // the sinc functions.
	        double sinc_x = 1.0, sinc_y=1.0, sinc_z=1.0;
	        double ax = M_PI*fx*dx;
	        double ay = M_PI*fy*dy;
	        double az = M_PI*fz*dz;
	        if (fx != 0.0) sinc_x = ax/sin(ax);
	        if (fy != 0.0) sinc_y = ay/sin(ay);
	        if (fz != 0.0) sinc_z = az/sin(az);
	        double grid_cor = pow(sinc_x*sinc_y*sinc_z, 2*interp_order);
          if (do_interlacing) grid_cor *= 0.25;

          double dkr = ddg[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
          double dki = ddg[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
          if (do_interlacing) {
            double kh = ax+ay+az;
            double sink = sin(kh);
            double cosk = cos(kh);
            double dkr_2 = ddg_2[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
            double dki_2 = ddg_2[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
            dkr += dkr_2*cosk - dki_2*sink;
            dki += dkr_2*sink + dki_2*cosk;
          }
 
          double prefac = 1.0;
          if (k == 0) {
            if (j == 0) {
              if ((i != 0) && (i != NX/2)) prefac = 0.5;
            } else {
              if (i != NX/2) prefac=0.5;
            } 
          }

          double L2 = 1.5*mu*mu - 0.5;
          double L4 = 4.375*mu*mu*mu*mu - 3.75*mu*mu + 0.375;
          double power = prefac*(dkr*dkr+dki*dki)*grid_cor;
          Pg0[kbin] += power;
          Pg2[kbin] += L2*power;
          Pg4[kbin] += L4*power;
          Nk[kbin]++;
          if (do_2D) {
            if ((bin >= 0) && (bin < NB*NMU)) {
              Pg0_2D[bin] += power;
              Nk_2D[bin]++;
            }
          }
	      }	
      }
    } 
  }
  
  // Sum the power and nmodes over the processors, then get processor 0 to do the
  // rest of the work
  double * Pg0_glob, * Pg2_glob, * Pg4_glob, * Nk_glob;
  if (ThisTask == 0) {
    Pg0_glob = (double*)calloc(NB, sizeof(double));       // binned power spectrum
    Pg2_glob = (double*)calloc(NB, sizeof(double));       // binned power spectrum
    Pg4_glob = (double*)calloc(NB, sizeof(double));       // binned power spectrum
    Nk_glob = (double*)calloc(NB, sizeof(double));       // number of modes in each bin
  }
  MPI_Reduce(Pg0, Pg0_glob, NB, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Pg2, Pg2_glob, NB, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Pg4, Pg4_glob, NB, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Nk, Nk_glob, NB, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (ThisTask == 0) {
    for(int i=0;i<NB;i++) {
      if(Nk[i]>0.0 && Pgfill[i]==0) {
        Pg0_glob[i]    -= Nk_glob[i]*nbwsq;
        Pg0_glob[i]    /= Nk_glob[i]*nbsqwsq; 
        Pg2_glob[i]    *= 5.0/(Nk_glob[i]*nbsqwsq);
        Pg4_glob[i]    *= 9.0/(Nk_glob[i]*nbsqwsq);
        Pgfill[i] = 1;
      }
    }
    
    // remove last bin -- this may only be part filled!
    for(int i=NB-1;i>=0;i--) {
      if(Pgfill[i]==1) {
        printf("last bin %d\n",i);
        Pgfill[i]=0; Pg0_glob[i]=0.0; Pg2_glob[i]=0.0; Pg4_glob[i]=0.0; Nk_glob[i]=0.0; break;
      }
    }   
   
    // output power spectrum values
    for(int i=0;i<NB;i++) {
      double kp;
      if (do_log) {
        kp = pow(10.0, log10(mink)+((float)i+0.5)*binwidth);
      } else {
        kp = mink+((float)i+0.5)*binwidth;
      }
      fprintf(fout,"%g %g %g %g %g\n",kp,Pg0_glob[i],Pg2_glob[i],Pg4_glob[i],Nk_glob[i]);
    }
    fclose(fout);
   
    free(Pg0_glob);
    free(Pg2_glob);
    free(Pg4_glob);
    free(Nk_glob);
  }

  if (do_2D) {
    // Sum the power and nmodes over the processors, then get processor 0 to do the
    // rest of the work
    double * Pg0_2D_glob, * Nk_2D_glob;
    if (ThisTask == 0) {
      Pg0_2D_glob = (double*)calloc(NB*NMU, sizeof(double));       // binned power spectrum
      Nk_2D_glob = (double*)calloc(NB*NMU, sizeof(double));       // number of modes in each bin
    }
    MPI_Reduce(Pg0_2D, Pg0_2D_glob, NB*NMU, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(Nk_2D, Nk_2D_glob, NB*NMU, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (ThisTask == 0) {
      for(int i=0;i<NB*NMU;i++) {
        if(Nk_2D[i]>0.0) {
          Pg0_2D_glob[i]    -= Nk_2D_glob[i]*nbwsq; 
          Pg0_2D_glob[i]    /= Nk_2D_glob[i]*nbsqwsq; 
        }
      }
     
      // output power spectrum values
      for(int i=0;i<NB;i++) {
        double kp;
        if (do_log) {
          kp = pow(10.0, log10(mink)+((float)i+0.5)*binwidth);
        } else {
          kp = mink+((float)i+0.5)*binwidth;
        }
        for(int j=0;j<NMU;j++) {
          double mup = ((float)j+0.5)/NMU;
          fprintf(fout_2D,"%g %g %g %g\n",kp,mup,Pg0_2D_glob[i*NMU+j],Nk_2D_glob[i*NMU+j]);
        }
      }
      fclose(fout_2D);
     
      free(Pg0_2D_glob);
      free(Nk_2D_glob);

    }
  }

  fftw_destroy_plan(dp_r2c);
  free(ddg);
  if (do_interlacing) {
    fftw_destroy_plan(dp_r2c_2);
    free(ddg_2);
  }
  free(Pg0);
  free(Pg2);
  free(Pg4);
  free(Pgfill);
  free(Nk);

  MPI_Finalize();

  return 0;
}

// This catches I/O errors occuring for fread(). In this case we better stop.
// ==========================================================================
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream) {
  size_t nread;
  if((nread = fread(ptr, size, nmemb, stream)) != nmemb) {
    printf("I/O error (fread) has occurred\n");
    FatalError((char *)"my_fread", 484);
  }
  return nread;
}

// Error message
// =============
void FatalError(char* filename, int linenum) {
  printf("Fatal Error at line %d in routine %s\n", linenum, filename);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

double add_to_grid(double x, double y, double z, double w, double xmin, double xmax, int nx, double dx, double dy, double dz, double * ddg) {

  if (interp_order == 1) {

    if (ThisTask == NTask-1) {
      if ((x < xmin) || (x >= xmax)) return 0.0;
    } else { 
      if ((x < xmin) || (x > xmax)) return 0.0;
    }

    int ix = (int)(((double)x-xmin)/dx);
    int iy = (int)(((double)y-YMIN)/dy);
    int iz = (int)(((double)z-ZMIN)/dz);
    if (ix == nx) ix = nx-1;
    if (iy == NY) iy = NY-1;
    if (iz == NZ) iz = NZ-1;

    unsigned long long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
    if (ind >= Total_size) {
      printf("%d, %lf, %lf, %lf, %llu\n", ThisTask, x, y, z, ind);
      FatalError("read_data", 133);
    }
    ddg[ind] += w;

  } else if (interp_order == 2) {

    if (ThisTask == NTask-1) {
      if ((x < xmin) || (x > xmax)) return 0.0;
    } else { 
      if ((x < xmin) || (x > xmax)) return 0.0;
    }

    double scalex = (x-xmin)/dx;
    double scaley = (y-YMIN)/dy;
    double scalez = (z-ZMIN)/dz;

    int ix = (int)scalex;
    int iy = (int)scaley;
    int iz = (int)scalez;
    double idx = scalex-(double)ix;
    double idy = scaley-(double)iy;
    double idz = scalez-(double)iz;
    double itx = 1.0-idx;
    double ity = 1.0-idy;
    double itz = 1.0-idz;
            
    int ixneigh = ix+1;
    int iyneigh = iy+1;
    int izneigh = iz+1;
    if(ixneigh >= nx) ixneigh = nx;
    if(iyneigh >= NY) {
      iyneigh = 0;
      if (iy >= NY) iy = NY-1;
    }
    if(izneigh >= NZ) {
      izneigh = 0;
      if (iz >= NZ) iz = NZ-1;
    }

    ddg[iz+2*(NZ/2+1)*(iy+NY*ix)]           += w*itx*ity*itz;
    ddg[izneigh+2*(NZ/2+1)*(iy+NY*ix)]      += w*itx*ity*idz;
    ddg[iz+2*(NZ/2+1)*(iyneigh+NY*ix)]      += w*itx*idy*itz;
    ddg[izneigh+2*(NZ/2+1)*(iyneigh+NY*ix)] += w*itx*idy*idz;

    ddg[iz+2*(NZ/2+1)*(iy+NY*ixneigh)]           += w*idx*ity*itz;
    ddg[izneigh+2*(NZ/2+1)*(iy+NY*ixneigh)]      += w*idx*ity*idz;
    ddg[iz+2*(NZ/2+1)*(iyneigh+NY*ixneigh)]      += w*idx*idy*itz;
    ddg[izneigh+2*(NZ/2+1)*(iyneigh+NY*ixneigh)] += w*idx*idy*idz;

  } else if (interp_order == 3) {

    if (x < dx) x += XMAX-XMIN;
    if (ThisTask == NTask-1) {
      if ((x < xmin+dx) || (x >= xmax+dx)) return 0.0;
    } else { 
      if ((x < xmin+dx) || (x > xmax+dx)) return 0.0;
    }

    double scalex = (x-xmin)/dx;
    double scaley = (y-YMIN)/dy;
    double scalez = (z-ZMIN)/dz;

    int ix = (int)floor(scalex+0.5);
    int iy = (int)floor(scaley+0.5);
    int iz = (int)floor(scalez+0.5);
    double idx = scalex-(double)ix;
    double idy = scaley-(double)iy;
    double idz = scalez-(double)iz;
    double itx = 0.75 - idx*idx;
    double ity = 0.75 - idy*idy;
    double itz = 0.75 - idz*idz;
    double isx = 0.5*(0.5-idx)*(0.5-idx);
    double isy = 0.5*(0.5-idy)*(0.5-idy);
    double isz = 0.5*(0.5-idz)*(0.5-idz);
    idx = 0.5*(0.5+idx)*(0.5+idx);
    idy = 0.5*(0.5+idy)*(0.5+idy);
    idz = 0.5*(0.5+idz)*(0.5+idz);

    int ixneighlow = ix-1;
    int iyneighlow = iy-1;
    int izneighlow = iz-1; 
    int ixneighhi = ix+1;
    int iyneighhi = iy+1;
    int izneighhi = iz+1; 
    //printf("%lf, %d, %d, %d, %lf\n", x, ixneighlow, ix, ixneighhi, dx);
    if (iyneighlow < 0) iyneighlow = NY-1;
    if (izneighlow < 0) izneighlow = NZ-1;
    if(ixneighhi >= nx) ixneighhi = nx;
    if(iyneighhi >= NY) {
      iyneighhi -= NY;
      if(iy >= NY) iy = 0;
      if(iyneighhi >= 1) iyneighhi = 1;
    }
    if(izneighhi >= NZ) {
      izneighhi -= NZ;
      if(iz >= NZ) iz = 0;
      if(izneighhi >= 1) izneighhi = 1;
    }

    ddg[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)] += w*isx*isy*isz;
    ddg[iz+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)]         += w*isx*isy*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)]  += w*isx*isy*idz;
    ddg[izneighlow+2*(NZ/2+1)*(iy+NY*ixneighlow)] += w*isx*ity*isz;
    ddg[iz+2*(NZ/2+1)*(iy+NY*ixneighlow)]         += w*isx*ity*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iy+NY*ixneighlow)]  += w*isx*ity*idz;
    ddg[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)] += w*isx*idy*isz;
    ddg[iz+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)]         += w*isx*idy*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)]  += w*isx*idy*idz;

    ddg[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ix)] += w*itx*isy*isz;
    ddg[iz+2*(NZ/2+1)*(iyneighlow+NY*ix)]         += w*itx*isy*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ix)]  += w*itx*isy*idz;
    ddg[izneighlow+2*(NZ/2+1)*(iy+NY*ix)] += w*itx*ity*isz;
    ddg[iz+2*(NZ/2+1)*(iy+NY*ix)]         += w*itx*ity*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iy+NY*ix)]  += w*itx*ity*idz;
    ddg[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ix)] += w*itx*idy*isz;
    ddg[iz+2*(NZ/2+1)*(iyneighhi+NY*ix)]         += w*itx*idy*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ix)]  += w*itx*idy*idz;

    ddg[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)] += w*idx*isy*isz;
    ddg[iz+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]         += w*idx*isy*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]  += w*idx*isy*idz;
    ddg[izneighlow+2*(NZ/2+1)*(iy+NY*ixneighhi)] += w*idx*ity*isz;
    ddg[iz+2*(NZ/2+1)*(iy+NY*ixneighhi)]         += w*idx*ity*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iy+NY*ixneighhi)]  += w*idx*ity*idz;
    ddg[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)] += w*idx*idy*isz;
    ddg[iz+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]         += w*idx*idy*itz;
    ddg[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]  += w*idx*idy*idz;
  }
  return w;
}
