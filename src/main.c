/* ==========================================================================*/
/*   Version 1.0.             Cullan Howlett,                                */
/*   Copyright (c) 2019       University of Queensland.                      */
/*                                                                           */
/*   This file is part of CAPOW.                                             */
/*                                                                           */
/*   CAPOW is free software: you can redistribute it and/or modify           */
/*   it under the terms of the GNU General Public License as published by    */
/*   the Free Software Foundation, either version 3 of the License, or       */
/*   (at your option) any later version.                                     */
/*                                                                           */
/*   CAPOW is distributed in the hope that it will be useful,                */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*   GNU General Public License for more details.                            */
/*                                                                           */
/*   You should have received a copy of the GNU General Public License       */
/*   along with CAPOW.  If not, see <http://www.gnu.org/licenses/>.          */
/* ==========================================================================*/

/* ======================================================*/
/* This file contains the main driver routine for CAPOW. */
/* ======================================================*/
     
#include "vars.h"
#include "proto.h"

int main(int argc, char **argv) {

  // Set up MPI
  // ==========
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  fftw_mpi_init();

  if(argc < 2) {
    if(ThisTask == 0) {
      fprintf(stdout, "Input parameters not found\n");
      fprintf(stdout, "Call with <ParameterFile>\n");
      fflush(stdout);
    }
    ierr = MPI_Finalize();
    exit(1);
  }

  // Read the run parameters and setup code
  // ======================================
  if (ThisTask == 0) {
    printf("\nReading Input Parameters and setting up CAPOW\n");
    printf("=============================================\n");
    fflush(stdout);
  }
  read_parameterfile(argv[1]);

  if (ThisTask == 0) {
    printf("\nAllocating grid\n"); 
    printf("===============\n");
    fflush(stdout);
  }
  create_grids();

  if (ThisTask == 0) {
    printf("\nReading input data and assigning to grid\n"); 
    printf("========================================");
    fflush(stdout);
  }
  double NGRID = read_data_serial_ascii();

  if (ThisTask == 0) {
    printf("\nComputing shot-noise, normalisation and preparing for FFT\n"); 
    printf("=========================================================\n");
    fflush(stdout);
  }

  // Sum over all processors to get the global quantities
  // needed for the mean density, shot-noise and normalisation
  double NGRID_TOT = 0;
  MPI_Allreduce(&NGRID, &NGRID_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double mean  = NGRID_TOT/((double)NX*(double)NY*(double)NZ);
  double nbar  = NGRID_TOT/((XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN));
  double shot  = NGRID_TOT;
  double norm  = nbar*NGRID_TOT;
  if (ThisTask == 0) {
    printf("mean          = %g\n",mean);
    printf("nbar          = %g\n",nbar);
    printf("shot-noise    = %g\n",shot);
    printf("normalisation = %g\n",norm);
    fflush(stdout);
  }

  // Subtract the mean density
  int i, j, k;
  #pragma omp parallel for private(i,j,k) collapse(3)
  for (i=0; i<Local_nx; i++) {
    for (j=0; j<NY; j++) {
      for (k=0; k<NZ; k++) {
        unsigned long long ind = k+2*(NZ/2+1)*(j+NY*i);
        ddg[ind] -= mean;
        if (DoInterlacing) ddg_2[ind] -= mean;
      }
    }
  }

  if (ThisTask == 0) {
    printf("\nFourier transforming overdensity field\n"); 
    printf("======================================\n");
    fflush(stdout);
  }
  fftw_execute(plan);
  if (DoInterlacing) fftw_execute(plan_2);

  if (ThisTask == 0) {
    printf("\nIterating over grid cells\n"); 
    printf("=========================\n");
    fflush(stdout);
  }
  compute_power();

  if (ThisTask == 0) {
    printf("\nOutputting power spectra\n"); 
    printf("========================\n");
    fflush(stdout);
  }
  output_power(shot, norm);

  if (ThisTask == 0) {
    printf("\nCleaning up\n"); 
    printf("===========\n");
    fflush(stdout);
  }
  MPI_Finalize();

  return 0;
}

// Compute the power spectrum from the Fourier transformed overdensity
// ===================================================================
void compute_power(void) {

  // Set up the bins for P(k) data 
  binwidth=(Maxk-Mink)/(double)NK; // bin width
  if (Output2D) { 
    Pk_2D     = (double*)calloc(NK*NMU, sizeof(double));
    Nmodes_2D = (int*)calloc(NK*NMU, sizeof(int));
  }
  Pk0 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Pk2 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Pk4 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Nmodes = (int*)calloc(NK, sizeof(int));        // number of modes in each bin
  
  // calculate Nyquist frequency
  double ny_x = (M_PI/dx);
  double ny_y = (M_PI/dy);
  double ny_z = (M_PI/dz);
  double min_nyquist=ny_x;
  if (ny_y<min_nyquist) min_nyquist=ny_y;
  if (ny_z<min_nyquist) min_nyquist=ny_z;

  // Loop over all cells on this processor.
  int i, j, k;
  double sx = 1.0/((float)NX*dx), sy = 1.0/((float)NY*dy), sz = 1.0/((float)NZ*dz);
  #pragma omp parallel for private(i,j,k) collapse(3) reduction(+:Pk0[:NK],Pk2[:NK],Pk4[:NK],Nmodes[:NK],Pk_2D[:NK*NMU],Nmodes_2D[:NK*NMU])
  for (i=Local_x_start; i<Local_x_start+Local_nx; i++) {
    for (j=0; j<NY; j++) {
      for (k=0; k<=NZ/2; k++) {
  
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
          if (LOS == 1) {
            mu = 2.0*M_PI*fx/fktot;
          } else if (LOS == 2) {
            mu = 2.0*M_PI*fy/fktot;
          } else {
            mu = 2.0*M_PI*fz/fktot;
          }
        } else {
          mu = 0.0;
        }
        int kbin;
        if (OutputLog) {
          kbin = (int)((float)(log10(fktot)-Mink)/binwidth);
        } else {
          kbin = (int)((float)(fktot-Mink)/binwidth);
        }
        int mubin = (int)(mu*(float)NMU);
        int bin = kbin*NMU + mubin;

        if ((kbin >= 0) && (kbin < NK) && (fktot <= min_nyquist)) {
    
          // set up correction for gridding - in effect we're
          // convolving the density field with a top-hat function in
          // each direction, so we're multiplying each Fourier mode by
          // a sinc function. To correct this, we therefore divide by
          // the sinc functions.
          double sinc_x = 1.0, sinc_y=1.0, sinc_z=1.0;
          double ax = M_PI*fx*dx;
          double ay = M_PI*fy*dy;
          double az = M_PI*fz*dz;
          if (fx != 0.0) sinc_x = sin(ax)/ax;
          if (fy != 0.0) sinc_y = sin(ay)/ay;
          if (fz != 0.0) sinc_z = sin(az)/az;
          double grid_cor = pow(1.0/(sinc_x*sinc_y*sinc_z), 2*InterpOrder);

          // Compute the real and imaginary parts of the overdensity field
          double dkr = ddg[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
          double dki = ddg[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
          if (DoInterlacing) {
            grid_cor *= 0.25;
            double kh = ax+ay+az;
            double sink = sin(kh);
            double cosk = cos(kh);
            double dkr_2 = ddg_2[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
            double dki_2 = ddg_2[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
            dkr += dkr_2*cosk - dki_2*sink;
            dki += dkr_2*sink + dki_2*cosk;
          }
 
          // This is necessary because we are only looping over half the grid in the z axis
          // so we need to fully account for the assumed hermitian symmetry in the data
          if (k == 0) {
            if (j == 0) {
              if ((i != 0) && (i != NX/2)) grid_cor *= 0.5;
            } else {
              if (i != NX/2) grid_cor *= 0.5;
            } 
          }

          // Get the power, multiply by the relevant Legendre polynomials and add to the bin
          double L2 = 1.5*mu*mu - 0.5;
          double L4 = 4.375*mu*mu*mu*mu - 3.75*mu*mu + 0.375;
          double power = (dkr*dkr+dki*dki)*grid_cor;
          Pk0[kbin]   += power;
          Pk2[kbin]   += L2*power;
          Pk4[kbin]   += L4*power;
          Nmodes[kbin]++;
          if (Output2D) {
            if ((bin >= 0) && (bin < NK*NMU)) {
              Pk_2D[bin]    += power;
              Nmodes_2D[bin]++;
            }
          }
        } 
      }
    } 
  }

  return;
}

// Output the power spectrum
// =========================
void output_power(double shot, double norm) {

  // Sum the power and number of modes over all processors.
  int * Nmodes_glob, * Nmodes_2D_glob;
  double * Pk0_glob, * Pk2_glob, * Pk4_glob, * Pk_2D_glob;
  if (ThisTask == 0) {
    Pk0_glob = (double*)calloc(NK, sizeof(double));
    Pk2_glob = (double*)calloc(NK, sizeof(double));
    Pk4_glob = (double*)calloc(NK, sizeof(double));
    Nmodes_glob = (int*)calloc(NK, sizeof(int));
    if (Output2D) {
      Pk_2D_glob = (double*)calloc(NK*NMU, sizeof(double));
      Nmodes_2D_glob = (int*)calloc(NK*NMU, sizeof(int));
    }
  }
  MPI_Reduce(Pk0, Pk0_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Pk2, Pk2_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Pk4, Pk4_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Nmodes, Nmodes_glob, NK, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (Output2D) {
    MPI_Reduce(Pk_2D, Pk_2D_glob, NK*NMU, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(Nmodes_2D, Nmodes_2D_glob, NK*NMU, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  // Now get processor 0 to write out the files
  if (ThisTask == 0) {

    FILE *fout, *fout_2D;
    char fout_name[2000], fout_name_2D[2000];
    sprintf(fout_name, "%s/%s.lpow", OutputDir, FileBase);
    if (ThisTask == 0) {
      if((fout=fopen(fout_name,"w"))==NULL) { printf("cannot open output file: %s\n", fout_name); FatalError("read_data", 142); }
      printf("Writing multipoles to file: %s\n",fout_name);
      fflush(stdout);
      if (Output2D) {
        sprintf(fout_name_2D, "%s_2D", fout_name);
        if((fout_2D=fopen(fout_name_2D,"w"))==NULL) { printf("cannot open output file: %s\n", fout_name_2D); FatalError("read_data", 142); }
        printf("Writing 2D power spectrum to file:%s\n",fout_name_2D);
        fflush(stdout);
      }
    }

    // Normalise and apply shot noise correction
    int * Pkfill = (int*)calloc(NK, sizeof(int));
    for(int i=0;i<NK;i++) {
      if (Nmodes_glob[i]>0.0 && Pkfill[i]==0) {
        Pk0_glob[i] -= Nmodes_glob[i]*shot; 
        Pk0_glob[i] /= Nmodes_glob[i]*norm; 
        Pk2_glob[i] *= 5.0/(Nmodes_glob[i]*norm);
        Pk4_glob[i] *= 9.0/(Nmodes_glob[i]*norm);
        Pkfill[i] = 1;
      }
      if (Output2D) {
        if(Nmodes_2D_glob[i]>0.0) {
          Pk_2D_glob[i] -= Nmodes_2D_glob[i]*shot; 
          Pk_2D_glob[i] /= Nmodes_2D_glob[i]*norm; 
        }
      }
    }
    
    // Remove the last bin -- this may only be partially filled!
    for(int i=NK-1;i>=0;i--) {
      if(Pkfill[i]==1) {
        double kmax;
        if (OutputLog) {
          kmax = pow(10.0, Mink+((float)i+0.5)*binwidth);
        } else {
          kmax = Mink+((float)i+0.5)*binwidth;
        }
        printf("Last full bin: %d (k=%g)\n",i, kmax);
        fflush(stdout);
        Pkfill[i]=0; Pk0_glob[i]=0.0; Pk2_glob[i]=0.0; Pk4_glob[i]=0.0; Nmodes_glob[i]=0.0; break;
      }
    }   
   
    // Output the power spectrum values
    // You may want to change this based on your choice of output format
    fprintf(fout, "# k, pk0, pk2, pk4, Nmodes\n");
    for(int i=0;i<NK;i++) {
      double kp;
      if (OutputLog) {
        kp = pow(10.0, Mink+((float)i+0.5)*binwidth);
      } else {
        kp = Mink+((float)i+0.5)*binwidth;
      }
      fprintf(fout,"%g %g %g %g %d\n",kp,Pk0_glob[i],Pk2_glob[i],Pk4_glob[i],Nmodes_glob[i]);
    }
    fclose(fout);
   
    if (Output2D) {
      fprintf(fout_2D, "# k, mu, pk, Nmodes\n");
      for(int i=0;i<NK;i++) {
        double kp;
        if (OutputLog) {
          kp = pow(10.0, Mink+((float)i+0.5)*binwidth);
        } else {
          kp = Mink+((float)i+0.5)*binwidth;
        }
        for(int j=0;j<NMU;j++) {
          double mup = ((float)j+0.5)/NMU;
          fprintf(fout_2D,"%g %g %g %d\n",kp,mup,Pk_2D_glob[i*NMU+j],Nmodes_2D_glob[i*NMU+j]);
        }
      }
      fclose(fout_2D);
    }

    free(Pk0_glob);
    free(Pk2_glob);
    free(Pk4_glob);
    free(Pkfill);
    free(Nmodes_glob);
    if (Output2D) {
      free(Pk_2D_glob);
      free(Nmodes_2D_glob);
    }
  }

  return;
}

// Create the grid(s) 
// ==================
void create_grids(void) {

  // Compute the grid spacing
  dx  = (XMAX-XMIN)/(double)NX;
  dy  = (YMAX-YMIN)/(double)NY;
  dz  = (ZMAX-ZMIN)/(double)NZ;

  // Compute some grid allocation variables
  ptrdiff_t alloc_local = fftw_mpi_local_size_3d(NX, NY, NZ/2+1, MPI_COMM_WORLD, &Local_nx, &Local_x_start);
  alloc_slice = 2*NY*(NZ/2+1);
  last_slice = Local_nx*alloc_slice;
  Local_nxtra = Local_nx+InterpOrder;
  Total_size = 2*alloc_local+InterpOrder*alloc_slice;

  // Set the neighbouring tasks
  Local_nx_table = (int *)malloc(sizeof(int) * NTask);
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

  // Allocate the grids and create the FFTW plan
  ddg = (double*)calloc(Total_size,sizeof(double));
  plan = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg,(fftw_complex*)ddg,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  if (DoInterlacing) {
    ddg_2 = (double*)calloc(Total_size,sizeof(double));
    plan_2 = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_2,(fftw_complex*)ddg_2,MPI_COMM_WORLD,FFTW_ESTIMATE);   
  }

  return;
}

// Destroy the grid(s)
// ===================
void destroy_grids(void) {
  fftw_destroy_plan(plan);
  free(ddg);
  if (DoInterlacing) {
    fftw_destroy_plan(plan_2);
    free(ddg_2);
  }
  free(Pk0);
  free(Pk2);
  free(Pk4);
  free(Nmodes);
  if (Output2D) {
    free(Pk_2D);
    free(Nmodes_2D);
  }
  return;
}

// Error message
// =============
void FatalError(char* filename, int linenum) {
  printf("Fatal Error at line %d in routine %s\n", linenum, filename);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}
