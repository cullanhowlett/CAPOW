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

  char datafile[2000], randfile[2000], outfile[2000];

  // Set up MPI
  // ==========
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  fftw_mpi_init();

  if(argc < 5) {
    if(ThisTask == 0) {
      fprintf(stdout, "Input parameters not found\n");
      fprintf(stdout, "Call with <ParameterFile, DataFile, RandomFile, OutputFile>\n");
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
  sprintf(datafile, "%s", argv[2]);
  sprintf(randfile, "%s", argv[3]);
  sprintf(outfile, "%s", argv[4]);

  if (ThisTask == 0) {
    printf("\nAllocating grid\n"); 
    printf("===============\n");
    fflush(stdout);
  }
  create_grids();

  if (ThisTask == 0) {
    printf("\nReading input data and assigning to grid\n"); 
    printf("========================================\n");
    fflush(stdout);
  }
  double NPERIODIC;
  unsigned long long NDATA, NRAND;
  if (Periodic) {
    NPERIODIC = read_periodic_serial_ascii(datafile);
  } else if (Survey) {
    REDMIN = 1.0e30, REDMAX = -1.0e30;
    if(!(data = (struct survey_data *) malloc(NOBJ_Max*sizeof(struct survey_data)))) { printf("Task %d unable to allocate memory for data\n", ThisTask);  FatalError("main", 72); }
    //sprintf(datafile, "%s/%s", InputDir, FileBase);
    NDATA = read_survey_serial_ascii(datafile, data, 0);
    //if (Momentum != 1) {
    if(!(randoms = (struct survey_data *) malloc(NOBJ_Max*sizeof(struct survey_data)))) { printf("Task %d unable to allocate memory for randoms\n", ThisTask);  FatalError("main", 72); }
    //sprintf(randfile, "%s/%s", InputDir, RandFileBase);
    NRAND = read_survey_serial_ascii(randfile, randoms, 1);
    //}
    if (NBAR_Column < 0) compute_nbar(0, NDATA, NRAND);
    if (FKP_Column < 0) {
      if (ThisTask == 0) printf("Computing FKP Weights...\n");
      compute_fkp(NDATA, data);
      //if (Momentum != 1) compute_fkp(NRAND, randoms);
      compute_fkp(NRAND, randoms);
    }
    assign_survey_data(NDATA, data, 1.0, 1);
  }

  if (ThisTask == 0) {
    printf("\nComputing shot-noise and normalisation\n"); 
    printf("======================================\n");
    fflush(stdout);
  }

  // Sum over all processors to get the global quantities
  // needed for the mean density, shot-noise and normalisation
  double alpha, shot, norm;
  if (Periodic) {
    double NPERIODIC_TOT = 0, NSQ_TOT = 0.0, VR_TOT = 0.0, VR2_TOT = 0.0, VR3_TOT = 0.0, VR4_TOT = 0.0;
    MPI_Allreduce(&NPERIODIC, &NPERIODIC_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (Momentum) {
   	  MPI_Allreduce(&nsq, &NSQ_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&vr_ave, &VR_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&vr2_ave, &VR2_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&vr3_ave, &VR3_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&vr4_ave, &VR4_TOT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    }
    alpha  = NPERIODIC_TOT/((double)NX*(double)NY*(double)NZ);
    double nbar  = NPERIODIC_TOT/((XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN));
    shot  = NPERIODIC_TOT;   //P00
    norm  = nbar*NPERIODIC_TOT;
    if (Momentum) {
      NSQ_TOT /= (double)NTask;
      VR_TOT /= (NSQ_TOT * (double)NTask);
      VR2_TOT /= (NSQ_TOT * (double)NTask);
      VR3_TOT /= (NSQ_TOT * (double)NTask);
      VR4_TOT /= (NSQ_TOT * (double)NTask);
    }
    if (ThisTask == 0) {
      printf("mean          = %g\n",alpha);
      printf("nbar          = %g\n",nbar);
      printf("shot-noise    = %g\n",shot);
      printf("normalisation = %g\n",norm);
      if (Momentum) {
      	printf("vr_ave = %g\n",VR_TOT);
        printf("vr2_ave = %g\n",VR2_TOT); 
        printf("vr3_ave = %g\n",VR3_TOT); 
        printf("vr4_ave = %g\n",VR4_TOT); 
      }
      fflush(stdout);
    }
    if (Momentum == 1) shot = VR_TOT * NPERIODIC_TOT;   //P01
    if (Momentum == 2) shot = VR2_TOT * NPERIODIC_TOT; //P11
    if (Momentum == 3) shot = VR2_TOT * NPERIODIC_TOT; //P02
    if (Momentum == 4) shot = VR3_TOT * NPERIODIC_TOT;  //P12
    if (Momentum == 5) shot = VR4_TOT * NPERIODIC_TOT;  //P22
    if (Momentum == 6) shot = VR3_TOT * NPERIODIC_TOT;  //P03
    if (Momentum == 7) shot = VR4_TOT * NPERIODIC_TOT;  //P13
    if (Momentum == 8) shot = VR4_TOT * NPERIODIC_TOT;  //P04
  } else if (Survey) {
    double data_nbw = 0.0, data_nbwsq = 0.0, data_nbsqwsq = 0.0, data_nbsqwsq_pv = 0.0, data_vr = 0.0, data_vrsq = 0.0; 
    double rand_nbw = 0.0, rand_nbwsq = 0.0, rand_nbsqwsq = 0.0, rand_nbsqwsq_pv = 0.0;
    for (unsigned long long i=0;i<NDATA;i++) {
      data_nbw     += data[i].weight;
      data_nbwsq   += data[i].weight*data[i].weight;
      data_nbsqwsq += data[i].weight*data[i].weight*data[i].nbar;
      if (Momentum) {
        data_vr += data[i].weight*data[i].weight_pv*data[i].pv;
        data_vrsq += data[i].weight*data[i].weight*data[i].pv*data[i].pv;
        data_nbsqwsq_pv += data[i].weight_pv*data[i].weight_pv*data[i].nbar;
      }
    }
    for (unsigned long long i=0;i<NRAND;i++) {
      rand_nbw     += randoms[i].weight;
      rand_nbwsq   += randoms[i].weight*randoms[i].weight;
      rand_nbsqwsq += randoms[i].weight*randoms[i].weight*randoms[i].nbar;
      if (Momentum == 2) {
        rand_nbsqwsq_pv += randoms[i].weight_pv*randoms[i].weight_pv*randoms[i].nbar;
      }
    }

    alpha = data_nbw/rand_nbw;
    rand_nbwsq   *= alpha*alpha;
    rand_nbsqwsq *= alpha;

    if (Momentum == 0) {
	  shot = data_nbwsq + rand_nbwsq;
	  norm = rand_nbsqwsq;
    } else if (Momentum == 1) {
      shot = data_vrsq;
      norm = rand_nbsqwsq;
    } else if (Momentum == 2) {
      shot = data_vr;
      norm = sqrt(rand_nbsqwsq*rand_nbsqwsq_pv);
    } else if (Momentum == 3) {
      shot = data_vrsq;
      norm = sqrt(rand_nbsqwsq*rand_nbsqwsq_pv);
    }

    if (ThisTask == 0) {
      if (Momentum == 0) {
	      printf("alpha = %g\n",alpha);
		    printf("shot-noise [data, randoms]    = %g, %g\n",data_nbwsq, rand_nbwsq);
		    printf("normalisation [data, randoms] = %g, %g\n",data_nbsqwsq, rand_nbsqwsq);
	    } else if (Momentum == 1) {
        printf("shot-noise [data]    = %g \n",data_vrsq);
        printf("normalisation [data, randoms] = %g, %g\n",data_nbsqwsq, rand_nbsqwsq);
	    } else {
	  	  printf("alpha = %g\n",alpha);
	  	  printf("shot-noise    = %g\n",data_vr);
		    printf("normalisation = %g\n",norm);
	    }
      fflush(stdout);
    }
  }

  if (ThisTask == 0) {
    printf("\nSubtracting off expected density\n"); 
    printf("================================\n");
    fflush(stdout);
  }
  if (Periodic) {
    // Subtract the mean density off of each cell
    // P00, P01, P11, P02, P12, P22, P03, P13, P04
    if ((Momentum == 0) | (Momentum == 1) | (Momentum == 3) | (Momentum == 6) | (Momentum == 8)) {
      int i, j, k;
      for (i=0; i<Local_nx; i++) {
        for (j=0; j<NY; j++) {
          for (k=0; k<NZ; k++) {
            unsigned long long ind = k+2*(NZ/2+1)*(j+NY*i);
            ddg[ind] -= alpha;
            if (DoInterlacing) ddg_interlace[ind] -= alpha;
          }
        }
      }
    }
  } else if (Survey) {
    // Subtract the random counts
    if (Momentum != 1) assign_survey_data(NRAND, randoms, -alpha, 0);

    // Copy across the extra slices from the task on the left and add it to the leftmost slices
    // of the task on the right. Skip over tasks without any slices.
    if (InterpOrder > 1) {
      double * temp_ddg = (double *)calloc(InterpOrder*alloc_slice, sizeof(double));
      ierr = MPI_Sendrecv(&(ddg[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                          &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
      for (int i=0;i<InterpOrder*alloc_slice;i++) ddg[i] += temp_ddg[i];
      if (DoInterlacing) {
        ierr = MPI_Sendrecv(&(ddg_interlace[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                            &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
        for (int i=0;i<InterpOrder*alloc_slice;i++) ddg_interlace[i] += temp_ddg[i];
      }
      free(temp_ddg);
    }
    if ((Momentum != 0) && (Momentum != 1)) {
      if (InterpOrder > 1) {
        double * temp_ddg = (double *)calloc(InterpOrder*alloc_slice, sizeof(double));
        ierr = MPI_Sendrecv(&(ddg_mom[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                            &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
        for (int i=0;i<InterpOrder*alloc_slice;i++) ddg_mom[i] += temp_ddg[i];
        if (DoInterlacing) {
          ierr = MPI_Sendrecv(&(ddg_mom_interlace[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                              &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
          for (int i=0;i<InterpOrder*alloc_slice;i++) ddg_mom_interlace[i] += temp_ddg[i];
        }
        free(temp_ddg);
      }
    }
  }

  if (Periodic) {
    if (ThisTask == 0) {
      printf("\nFourier transforming overdensity field\n"); 
      printf("======================================\n");
      fflush(stdout);
    }
    fftw_execute(plan);
    if (DoInterlacing) fftw_execute(plan_interlace);
    if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)){
   	  fftw_execute(plan_mom);
      if (DoInterlacing) fftw_execute(plan_mom_interlace);
    }

    if (ThisTask == 0) {
      printf("\nIterating over grid cells\n"); 
      printf("=========================\n");
      fflush(stdout);
    }
    compute_periodic_power();
  } else if (Survey) {
    compute_survey_power();
  }

  if (ThisTask == 0) {
    printf("\nOutputting power spectra\n"); 
    printf("========================\n");
    fflush(stdout);
  }
  output_power(shot, norm, outfile);

  if (ThisTask == 0) {
    printf("\nCleaning up\n"); 
    printf("===========\n");
    fflush(stdout);
  }
  MPI_Finalize();

  return 0;
}

// Compute the power spectrum from the Fourier transformed overdensity for a periodic simulation
// =============================================================================================
void compute_periodic_power(void) {

  // Set up the bins for P(k) data 
  binwidth=(Maxk-Mink)/(double)NK; // bin width
  if (Output2D) { 
    Pk_2D     = (double*)calloc(NK*NMU, sizeof(double));
    Nmodes_2D = (int*)calloc(NK*NMU, sizeof(int));
  }
  Pk0 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Pk2 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Pk4 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  if (Odd_Multipoles) {
    Pk1 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
    Pk3 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  }
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
  //#pragma omp parallel for private(i,j,k) collapse(3) reduction(+:Pk0[:NK],Pk2[:NK],Pk4[:NK],Nmodes[:NK],Pk_2D[:NK*NMU],Nmodes_2D[:NK*NMU])
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
            double dkr_interlace = ddg_interlace[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
            double dki_interlace = ddg_interlace[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
            dkr += dkr_interlace*cosk - dki_interlace*sink;
            dki += dkr_interlace*sink + dki_interlace*cosk;
          }
          double dkr_mom, dki_mom;
          // P00, P01, P11, P02, P12, P22, P03, P13, P04
          if ((Momentum == 1) | (Momentum == 3) | (Momentum == 4) | (Momentum >= 6)) {
          	dkr_mom = ddg_mom[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
            dki_mom = ddg_mom[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
            if (DoInterlacing) {
              double kh = ax+ay+az;
              double sink = sin(kh);
              double cosk = cos(kh);
              double dkr_mom_interlace = ddg_mom_interlace[(2*k  )+2*(NZ/2+1)*(j+NY*(i-Local_x_start))];
              double dki_mom_interlace = ddg_mom_interlace[(2*k+1)+2*(NZ/2+1)*(j+NY*(i-Local_x_start))]; 
              dkr_mom += dkr_mom_interlace*cosk - dki_mom_interlace*sink;
              dki_mom += dkr_mom_interlace*sink + dki_mom_interlace*cosk;
            }
          }
 
          // This is necessary because we are only looping over half the grid in the z axis
          // so we need to fully account for the assumed hermitian symmetry in the data
          int NM = 2;
          if (k == 0) {
            if (j == 0) {
              if ((i != 0) && (i != NX/2)) NM = 1;
            } else {
              if (i != NX/2) NM = 1;
            } 
          }

          // Get the power, multiply by the relevant Legendre polynomials and add to the bin
          double power;
          double L2 = 1.5*mu*mu - 0.5;
          double L4 = 4.375*mu*mu*mu*mu - 3.75*mu*mu + 0.375;
          if ((Momentum == 0) | (Momentum == 2) | (Momentum == 5)) {
            power = (dkr*dkr+dki*dki)*grid_cor;
       	  } else if ((Momentum == 1) | (Momentum == 4) | (Momentum == 6)) {
       	  	power = (dkr*dki_mom - dki*dkr_mom)*grid_cor;
       	  } else {
            power = (dkr*dkr_mom + dki*dki_mom)*grid_cor;
          }
          Pk0[kbin] += NM*power;
          Pk2[kbin] += NM*L2*power;
          Pk4[kbin] += NM*L4*power;
          if (Odd_Multipoles) {
            double L3 = 2.5*mu*mu*mu - 1.5*mu;
            Pk1[kbin] += NM*mu*power;
            Pk3[kbin] += NM*L3*power;
          }


          Nmodes[kbin] += NM;
          if (Output2D) {
            int mubin = (int)(((mu+1.0)/2.0)*(float)NMU);
            if (mubin >= NMU) mubin = NMU-1;
            int bin = kbin*NMU + mubin;
            if ((bin >= 0) && (bin < NK*NMU)) {
              Pk_2D[bin]    += NM*power;
              Nmodes_2D[bin] += NM;
            }
          }
        } 
      }
    } 
  }

  return;
}

// Compute the power spectrum from the Fourier transformed overdensity for a survey
// ===============================================================================
void compute_survey_power(void) {

  // Set up the bins for P(k) data 
  binwidth=(Maxk-Mink)/(double)NK; // bin width
  if (Output2D) { 
    Pk_2D     = (double*)calloc(NK*NMU, sizeof(double));
    Nmodes_2D = (int*)calloc(NK*NMU, sizeof(int));
  }
  Pk0 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Pk2 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  Pk4 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  if (Odd_Multipoles) {
    Pk1 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
    Pk3 = (double*)calloc(NK, sizeof(double));       // binned power spectrum
  }
  Nmodes = (int*)calloc(NK, sizeof(int));        // number of modes in each bin

  // Loop over the multipoles
  int multipole_counter = 2;
  if (Odd_Multipoles) multipole_counter = 1;
  for (int multipole=0; multipole<=4; multipole+=multipole_counter) {

    // For the monopole, we only need to take the Fourier transform
    // and loop over the grid
    if (multipole == 0) {
      if (ThisTask == 0) {
        printf("\nFourier transforming monopole overdensity field\n"); 
        printf("===============================================\n");
        fflush(stdout);
      }
      
      // We need to do this as we store the density field before the FFT
	  for (int ix=0; ix<Local_nx; ix++) {
	    for (int iy=0; iy<NY; iy++) {
	      for (int iz=0; iz<NZ; iz++) {
	        long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
	        ddg_2[ind] = ddg[ind];
	        if (DoInterlacing) ddg_interlace_2[ind] = ddg_interlace[ind];
	      }
	    }
	  }

      fftw_execute(plan);
      if (DoInterlacing) fftw_execute(plan_interlace);
      if ((Momentum != 0) && (Momentum != 1)) {
        for (int ix=0; ix<Local_nx; ix++) {
	      for (int iy=0; iy<NY; iy++) {
	        for (int iz=0; iz<NZ; iz++) {
	          long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
	          ddg_mom_2[ind] = ddg_mom[ind];
	          if (DoInterlacing) ddg_mom_interlace_2[ind] = ddg_mom_interlace[ind];
	        }
	      }
	    }

        fftw_execute(plan_mom);
        if (DoInterlacing) fftw_execute(plan_mom_interlace);
      }
      assign_survey_power(0, 0, 0, 0);
    }

    // For the higher order multipoles, we need to loop over directions of the 
    // unit vectors, multiply the density grid, do the Fourier transform and loop over grid cells

    // Dipole
    if (multipole == 1) {
      if (ThisTask == 0) {
        printf("\nFourier transforming dipole overdensity field\n"); 
        printf("=================================================\n");
        fflush(stdout);
      }

      for (int ii=0; ii<3; ii++) {
        if (ThisTask == 0) {
          printf("Indices: %d\n", ii); 
          fflush(stdout);
        }

        // Multiply the density grid by the real space unit vectors
        for (int ix=0; ix<Local_nx; ix++) {
          for (int iy=0; iy<NY; iy++) {
            for (int iz=0; iz<NZ; iz++) {
              long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
              double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
              double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
              if (r2 == 0.0) r2 = 1.0;
              ddg_2[ind] = ddg[ind]*vec[ii]/sqrt(r2);
              if (DoInterlacing) ddg_interlace_2[ind] = ddg_interlace[ind]*vec[ii]/sqrt(r2);
            }
          }
        }

        fftw_execute(plan);
        if (DoInterlacing) fftw_execute(plan_interlace);
        if ((Momentum != 0) && (Momentum != 1)) {
          for (int ix=0; ix<Local_nx; ix++) {
	        for (int iy=0; iy<NY; iy++) {
	          for (int iz=0; iz<NZ; iz++) {
                long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
                double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                if (r2 == 0.0) r2 = 1.0;
                ddg_mom_2[ind] = ddg_mom[ind]*vec[ii]/sqrt(r2);
                if (DoInterlacing) ddg_mom_interlace_2[ind] = ddg_mom_interlace[ind]*vec[ii]/sqrt(r2);
	          }
	        }
	      }

          fftw_execute(plan_mom);
          if (DoInterlacing) fftw_execute(plan_mom_interlace);
        }
        assign_survey_power(1, ii, 0, 0);
      }
    }

    // Quadrupole
    if (multipole == 2) {
      if (ThisTask == 0) {
        printf("\nFourier transforming quadrupole overdensity field\n"); 
        printf("=================================================\n");
        fflush(stdout);
      }

      for (int ii=0; ii<3; ii++) {
        for (int jj=ii; jj<3; jj++) {
          if (ThisTask == 0) {
            printf("Indices: %d %d\n", ii, jj); 
            fflush(stdout);
          }

          // Multiply the density grid by the real space unit vectors
          for (int ix=0; ix<Local_nx; ix++) {
            for (int iy=0; iy<NY; iy++) {
              for (int iz=0; iz<NZ; iz++) {
                long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
                double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                if (r2 == 0.0) r2 = 1.0;
                ddg_2[ind] = ddg[ind]*vec[ii]*vec[jj]/r2;
                if (DoInterlacing) ddg_interlace_2[ind] = ddg_interlace[ind]*vec[ii]*vec[jj]/r2;
              }
            }
          }

          fftw_execute(plan);
          if (DoInterlacing) fftw_execute(plan_interlace);

          if ((Momentum != 0) && (Momentum != 1)) {
            for (int ix=0; ix<Local_nx; ix++) {
	          for (int iy=0; iy<NY; iy++) {
	            for (int iz=0; iz<NZ; iz++) {
                  long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                  double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
                  double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                  if (r2 == 0.0) r2 = 1.0;
                  ddg_mom_2[ind] = ddg_mom[ind]*vec[ii]*vec[jj]/r2;
                  if (DoInterlacing) ddg_mom_interlace_2[ind] = ddg_mom_interlace[ind]*vec[ii]*vec[jj]/r2;
	            }
	          }
	        }

            fftw_execute(plan_mom);
            if (DoInterlacing) fftw_execute(plan_mom_interlace);
          }
          assign_survey_power(2, ii, jj, 0);
        }
      }
    }

    // Octopole
    if (multipole == 3) {
      if (ThisTask == 0) {
        printf("\nFourier transforming octopole overdensity field\n"); 
        printf("=================================================\n");
        fflush(stdout);
      }

      for (int ii=0; ii<3; ii++) {
        int jjlow=ii, jjhigh=ii+1;
        if (ii == 0) jjhigh = ii+2;
        for (int jj=jjlow; jj<jjhigh; jj++) {
          int kklow=0, kkhigh=3;
          if (jj != ii) kklow = 2;
          for (int kk=kklow; kk<kkhigh; kk++) {
            if (ThisTask == 0) {
              printf("Indices: %d %d %d\n", ii, jj, kk); 
              fflush(stdout);
            }

            // Multiply the density grid by the real space unit vectors
            for (int ix=0; ix<Local_nx; ix++) {
              for (int iy=0; iy<NY; iy++) {
                for (int iz=0; iz<NZ; iz++) {
                  long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                  double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
                  double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                  if (r2 == 0.0) r2 = 1.0;
                  ddg_2[ind] = ddg[ind]*vec[ii]*vec[jj]*vec[kk]/(r2*sqrt(r2));
                  if (DoInterlacing) ddg_interlace_2[ind] = ddg_interlace[ind]*vec[ii]*vec[jj]*vec[kk]/(r2*sqrt(r2));
                }
              }
            }

            fftw_execute(plan);
            if (DoInterlacing) fftw_execute(plan_interlace);

            if ((Momentum != 0) && (Momentum != 1)) {
              for (int ix=0; ix<Local_nx; ix++) {
	            for (int iy=0; iy<NY; iy++) {
	              for (int iz=0; iz<NZ; iz++) {
                    long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                    double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
                    double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                    if (r2 == 0.0) r2 = 1.0;
                    ddg_mom_2[ind] = ddg_mom[ind]*vec[ii]*vec[jj]*vec[kk]/(r2*sqrt(r2));
                    if (DoInterlacing) ddg_mom_interlace_2[ind] = ddg_mom_interlace[ind]*vec[ii]*vec[jj]*vec[kk]/(r2*sqrt(r2));
	              }
	            }
	          }

              fftw_execute(plan_mom);
              if (DoInterlacing) fftw_execute(plan_mom_interlace);
            }
            assign_survey_power(3, ii, jj, kk);
          }
        }
      }
    }

    // Hexadecapole
    if (multipole == 4) {
      if (ThisTask == 0) {
        printf("\nFourier transforming hexadecapole overdensity field\n"); 
        printf("===================================================\n");
        fflush(stdout);
      }


      for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
          int kklow=0, kkhigh=3;
          if (jj != ii) {
            if (ii == 1) {
              kklow = 2;
            } else if (ii == 2) {
              if (jj == 1) {
                continue;
              } else {
                kklow = 1;
                kkhigh = 2;
              } 
            } else {
              kklow = jj;
            }
          }
          for (int kk=kklow; kk<kkhigh; kk++) {
            if (ThisTask == 0) {
              printf("Indices: %d %d %d\n", ii, jj, kk); 
              fflush(stdout);
            }

            // Multiply the density grid by the real space unit vectors
            for (int ix=0; ix<Local_nx; ix++) {
              for (int iy=0; iy<NY; iy++) {
                for (int iz=0; iz<NZ; iz++) {
                  long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                  double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, (iy)*dy+YMIN-Y_Origin, (iz)*dz+ZMIN-Z_Origin};
                  double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                  if (r2 == 0.0) r2 = 1.0;
                  ddg_2[ind] = ddg[ind]*vec[ii]*vec[ii]*vec[jj]*vec[kk]/(r2*r2);
                  if (DoInterlacing) ddg_interlace_2[ind] = ddg_interlace[ind]*vec[ii]*vec[ii]*vec[jj]*vec[kk]/(r2*r2);
                }
              }
            }

            fftw_execute(plan);
            if (DoInterlacing) fftw_execute(plan_interlace);

            if ((Momentum != 0) && (Momentum != 1)) {
              for (int ix=0; ix<Local_nx; ix++) {
	            for (int iy=0; iy<NY; iy++) {
	              for (int iz=0; iz<NZ; iz++) {
                    long ind = iz+2*(NZ/2+1)*(iy+NY*ix);
                    double vec[3] = {(ix+Local_x_start)*dx+XMIN-X_Origin, iy*dy+YMIN-Y_Origin, iz*dz+ZMIN-Z_Origin};
                    double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                    if (r2 == 0.0) r2 = 1.0;
                    ddg_mom_2[ind] = ddg_mom[ind]*vec[ii]*vec[ii]*vec[jj]*vec[kk]/(r2*r2);
                    if (DoInterlacing) ddg_mom_interlace_2[ind] = ddg_mom_interlace[ind]*vec[ii]*vec[ii]*vec[jj]*vec[kk]/(r2*r2);
	              }
	            }
	          }

              fftw_execute(plan_mom);
              if (DoInterlacing) fftw_execute(plan_mom_interlace);
            }

            assign_survey_power(4, ii, jj, kk);
          }
        }
      }
    } 
  }

  return;
}

// Loop over the grid, compute the correct power and assign to the
// relevant multipole bin based on the multipole we are interested in.
// ====================================================================
void assign_survey_power(int multipole, int ii, int jj, int kk) {

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
        int kbin;
        if (OutputLog) {
          kbin = (int)((float)(log10(fktot)-Mink)/binwidth);
        } else {
          kbin = (int)((float)(fktot-Mink)/binwidth);
        }

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

          // Compute the power
          double power;
          double kvec[3] = {2.*M_PI*fx, 2.*M_PI*fy, 2.*M_PI*fz};
          long long ind = (2*k)+2*(NZ/2+1)*(j+NY*(i-Local_x_start));
          double dkr = ddg_2[ind];
          double dki = ddg_2[ind+1]; 
          if (DoInterlacing) {
            grid_cor *= 0.25;
            double kh = ax+ay+az;
            double sink = sin(kh);
            double cosk = cos(kh);
            double dkr_interlace = ddg_interlace_2[ind];
            double dki_interlace = ddg_interlace_2[ind+1]; 
            dkr += dkr_interlace*cosk - dki_interlace*sink;
            dki += dkr_interlace*sink + dki_interlace*cosk;
          }
          double dkr_mom, dki_mom;
          if ((Momentum != 0) && (Momentum != 1)) {
          	dkr_mom = ddg_mom_2[ind];
            dki_mom = ddg_mom_2[ind+1]; 
            if (DoInterlacing) {
              double kh = ax+ay+az;
              double sink = sin(kh);
              double cosk = cos(kh);
              double dkr_mom_interlace_2 = ddg_mom_interlace_2[ind];
              double dki_mom_interlace_2 = ddg_mom_interlace_2[ind+1]; 
              dkr_mom += dkr_mom_interlace_2*cosk - dki_mom_interlace_2*sink;
              dki_mom += dkr_mom_interlace_2*sink + dki_mom_interlace_2*cosk;
            }
          }
          if ((Momentum == 0) || (Momentum == 1)) {
	        if (multipole == 0) {
	          F0[ind] = dkr;
	          F0[ind+1] = dki;
	          power = (dkr*dkr+dki*dki)*grid_cor;
	        } else if (multipole % 2) {
	          power = (dkr*F0[ind+1]-dki*F0[ind])*grid_cor;
	        } else {
	          power = (dkr*F0[ind]+dki*F0[ind+1])*grid_cor;
	        }
	      } else {
	        if (multipole == 0) {
	          F0[ind] = dkr;
	          F0[ind+1] = dki;
	          F0_mom[ind] = dkr_mom;
	          F0_mom[ind+1] = dki_mom;
	          power = (dkr*dki_mom-dki*dkr_mom)*grid_cor;
	        } else if (multipole % 2) {
	          power = 0.5*(dkr*F0_mom[ind] + dki*F0_mom[ind+1] + dkr_mom*F0[ind] + dki_mom*F0[ind+1])*grid_cor;
	        } else {
	          power = 0.5*(dkr*F0_mom[ind+1] - dki*F0_mom[ind] + dkr_mom*F0[ind+1] - dki_mom*F0[ind])*grid_cor;
	        }
	      }

          // This is necessary because we are only looping over half the grid in the z axis
          // so we need to fully account for the assumed hermitian symmetry in the data
          int NM = 2;
          if (k == 0) {
            if (j == 0) {
              if ((i != 0) && (i != NX/2)) NM = 1;
            } else {
              if (i != NX/2) NM = 1;
            } 
          }

          // Add the power to the bin
          power *= NM;
          if (multipole == 0) {
            Pk0[kbin] += power;
            Pk2[kbin] -= 0.5*power;
            Pk4[kbin] += 0.375*power;
            Nmodes[kbin] += NM;
          }
          if (multipole == 1) {
            if (fktot > 0.0) {
              power *= kvec[ii]/fktot;
              Pk1[kbin] += power;
              Pk3[kbin] -= 1.5*power;
            }
          }
          if (multipole == 2) {
            if (fktot > 0.0) {
              double kprefac = 1.0;
              if (ii != jj) {
                kprefac = 2.0;
              }
              power *= kprefac*kvec[ii]*kvec[jj]/(fktot*fktot);
              Pk2[kbin] += 1.5*power;
              Pk4[kbin] -= 3.75*power;
            }
          }
          if (multipole == 3) {
            if (fktot > 0.0) {
              double kprefac = 1.0;
              if (ii != kk) {
                if (ii != jj) {
                  kprefac = 6.0;
                } else {
                  kprefac = 3.0;
                }
              }
              Pk3[kbin] += 2.5*kprefac*kvec[ii]*kvec[jj]*kvec[kk]/(fktot*fktot*fktot)*power;
            }
          }
          if (multipole == 4) {
            if (fktot > 0) {
              double kprefac = 1.0;
              if (ii == jj) {
                if (ii != kk) kprefac = 4.0;
              } else {
                if (jj == kk) {
                  kprefac = 6.0;
                } else {
                  kprefac = 12.0;
                }
              }
              Pk4[kbin] += 4.375*kprefac*kvec[ii]*kvec[ii]*kvec[jj]*kvec[kk]/(fktot*fktot*fktot*fktot)*power;
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
void output_power(double shot, double norm, char* fout_name) {

  // Sum the power and number of modes over all processors.
  int * Nmodes_glob, * Nmodes_2D_glob;
  double * Pk0_glob, * Pk1_glob, * Pk2_glob, * Pk3_glob, * Pk4_glob, * Pk_2D_glob;
  if (ThisTask == 0) {
    Pk0_glob = (double*)calloc(NK, sizeof(double));
    Pk2_glob = (double*)calloc(NK, sizeof(double));
    Pk4_glob = (double*)calloc(NK, sizeof(double));
    Nmodes_glob = (int*)calloc(NK, sizeof(int));
    if (Odd_Multipoles) {
      Pk1_glob = (double*)calloc(NK, sizeof(double));
      Pk3_glob = (double*)calloc(NK, sizeof(double));
    }
    if (Output2D) {
      Pk_2D_glob = (double*)calloc(NK*NMU, sizeof(double));
      Nmodes_2D_glob = (int*)calloc(NK*NMU, sizeof(int));
    }
  }

  MPI_Reduce(Pk0, Pk0_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Pk2, Pk2_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Pk4, Pk4_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Nmodes, Nmodes_glob, NK, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (Odd_Multipoles) {
    MPI_Reduce(Pk1, Pk1_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(Pk3, Pk3_glob, NK, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  }
  if (Output2D) {
    MPI_Reduce(Pk_2D, Pk_2D_glob, NK*NMU, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(Nmodes_2D, Nmodes_2D_glob, NK*NMU, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  // Now get processor 0 to write out the files
  if (ThisTask == 0) {

    FILE *fout, *fout_2D;
    char fout_name_2D[2000];
  	if((fout=fopen(fout_name,"w"))==NULL) { printf("cannot open output file: %s\n", fout_name); FatalError("read_data", 142); }
  	printf("Writing multipoles to file: %s\n",fout_name);
  	fflush(stdout);
  	if (Output2D) {
  	  sprintf(fout_name_2D, "%s_2D", fout_name);
  	  if((fout_2D=fopen(fout_name_2D,"w"))==NULL) { printf("cannot open output file: %s\n", fout_name_2D); FatalError("read_data", 142); }
  	  printf("Writing 2D power spectrum to file:%s\n",fout_name_2D);
  	  fflush(stdout);
  	}

    // Normalise and apply shot noise correction
    int * Pkfill = (int*)calloc(NK, sizeof(int));
    for(int i=0;i<NK;i++) {
      if (Nmodes_glob[i]>0.0 && Pkfill[i]==0) {
        Pk0_glob[i] -= Nmodes_glob[i]*shot; 
        Pk0_glob[i] /= Nmodes_glob[i]*norm; 
        Pk2_glob[i] *= 5.0/(Nmodes_glob[i]*norm);
        Pk4_glob[i] *= 9.0/(Nmodes_glob[i]*norm);
        if ((Odd_Multipoles)) {
          Pk1_glob[i] *= 3.0/(Nmodes_glob[i]*norm);
          Pk3_glob[i] *= 7.0/(Nmodes_glob[i]*norm);
        }
        Pkfill[i] = 1;
      }
    }
    if (Output2D) {
      for(int i=0;i<NK*NMU;i++) {
        if(Nmodes_2D_glob[i]>0.0) {
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
        if (Odd_Multipoles) { Pk1_glob[i]=0.0; Pk3_glob[i]=0.0; }
        Pkfill[i]=0; Pk0_glob[i]=0.0; Pk2_glob[i]=0.0; Pk4_glob[i]=0.0; Nmodes_glob[i]=0.0; break;
      }
    }   
   
    // Output the power spectrum values
    // You may want to change this based on your choice of output format
    if (Odd_Multipoles) {
      fprintf(fout, "# shot = %g\n", shot/norm);
      fprintf(fout, "# norm = %g\n", norm);
      fprintf(fout, "# k  pk0  pk1  pk2  pk3  pk4  Nmodes\n");
    } else {
      fprintf(fout, "# shot = %g\n", shot/norm);
      fprintf(fout, "# norm = %g\n", norm);
      fprintf(fout, "# k  pk0  pk2  pk4  Nmodes\n");
    }
    for(int i=0;i<NK;i++) {
      double kp;
      if (OutputLog) {
        kp = pow(10.0, Mink+((float)i+0.5)*binwidth);
      } else {
        kp = Mink+((float)i+0.5)*binwidth;
      }
      if (Odd_Multipoles) {
        fprintf(fout,"%g %g %g %g %g %g %d\n",kp,Pk0_glob[i],Pk1_glob[i],Pk2_glob[i],Pk3_glob[i],Pk4_glob[i],Nmodes_glob[i]);
      } else {
        fprintf(fout,"%g %g %g %g %d\n",kp,Pk0_glob[i],Pk2_glob[i],Pk4_glob[i],Nmodes_glob[i]);
      }
    }
    fclose(fout);
   
    if (Output2D) {
      fprintf(fout_2D, "# shot = %g\n", shot/norm);
      fprintf(fout_2D, "# norm = %g\n", norm);
      fprintf(fout_2D, "# k, mu, pk, Nmodes\n");
      for(int i=0;i<NK;i++) {
        double kp;
        if (OutputLog) {
          kp = pow(10.0, Mink+((float)i+0.5)*binwidth);
        } else {
          kp = Mink+((float)i+0.5)*binwidth;
        }
        for(int j=0;j<NMU;j++) {
          double mup = 2.0 * ((float)j+0.5)/NMU - 1.0;
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
    if (Odd_Multipoles) {
      free(Pk1_glob);
      free(Pk3_glob);
    }
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
  //printf("%ld, %ld, %ld\n", alloc_slice, last_slice, Total_size);

  if (ThisTask == 0) {
    printf("Setting neighbours...\n"); 
    fflush(stdout);
  }

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
      if (Periodic) {
        if (LeftTask < 0) {
          if (LOS == 0) {
            LeftTask = NTask - 1;
          } else {
            LeftTask = MPI_PROC_NULL;
            break;
          }
        }
      } else if (Survey) {
        if(LeftTask < 0) {
          LeftTask = MPI_PROC_NULL;
          break;
        }
      }
    } while(Local_nx_table[LeftTask] == 0);
      
    RightTask = ThisTask;
    do {
      RightTask++;
      if (Periodic) {
        if(RightTask >= NTask) {
          if (LOS == 0) {
            RightTask = 0;
          } else {
            RightTask = MPI_PROC_NULL;
            break;
          }
        }
      } else if (Survey) {
        if(RightTask >= NTask) {
          RightTask = MPI_PROC_NULL;
          break;
        }
      }
    } while(Local_nx_table[RightTask] == 0);
  }

  if (ThisTask == 0) {
    printf("Reserving memory...\n"); 
    fflush(stdout);
  }

  // Allocate the grids and create the FFTW plan
  ddg = (double*)calloc(Total_size,sizeof(double));
  if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) ddg_mom = (double*)calloc(Total_size,sizeof(double));
  if (Periodic) {
    plan = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg,(fftw_complex*)ddg,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) plan_mom = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_mom,(fftw_complex*)ddg_mom,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  } else if (Survey) {
    ddg_2 = (double*)calloc(Total_size,sizeof(double));
    plan = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_2,(fftw_complex*)ddg_2,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    F0 = (double*)malloc(Total_size*sizeof(double));
    if ((Momentum != 0) && (Momentum != 1)) {
      ddg_mom_2 = (double*)calloc(Total_size,sizeof(double));
      plan_mom = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_mom_2,(fftw_complex*)ddg_mom_2,MPI_COMM_WORLD,FFTW_ESTIMATE); 
      F0_mom = (double*)malloc(Total_size*sizeof(double));
    }
  }
  if (DoInterlacing) {
    ddg_interlace = (double*)calloc(Total_size,sizeof(double));
    if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) ddg_mom_interlace = (double*)calloc(Total_size,sizeof(double));
    if (Periodic) {
      plan_interlace = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_interlace,(fftw_complex*)ddg_interlace,MPI_COMM_WORLD,FFTW_ESTIMATE);   
      if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) plan_mom_interlace = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_mom_interlace,(fftw_complex*)ddg_mom_interlace,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    } else if (Survey) {
      ddg_interlace_2 = (double*)calloc(Total_size,sizeof(double));
      plan_interlace = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_interlace_2,(fftw_complex*)ddg_interlace_2,MPI_COMM_WORLD,FFTW_ESTIMATE);   
      if ((Momentum != 0) && (Momentum != 1)) {
    	ddg_mom_interlace_2 = (double*)calloc(Total_size,sizeof(double));
    	plan_mom_interlace = fftw_mpi_plan_dft_r2c_3d(NX,NY,NZ,ddg_mom_interlace_2,(fftw_complex*)ddg_mom_interlace_2,MPI_COMM_WORLD,FFTW_ESTIMATE); 
      }
    }
  }

  if (ThisTask == 0) {
    printf("Creating redshift-distance lookup table...\n"); 
    fflush(stdout);
  }

  // Generate a redshift-distance lookup table if necessary
  int nbins = 10000;
  REDMIN =  0.0;
  REDMAX = 10.0;
  double redbinwidth = (REDMAX-REDMIN)/(double)(nbins-1);
  double * ztemp = (double *)malloc(nbins*sizeof(double));
  double * rtemp = (double *)malloc(nbins*sizeof(double));
  ztemp[0] = 0.0, rtemp[0] = 0.0;
  for (int i=1;i<nbins;i++) {
      ztemp[i] = i*redbinwidth+REDMIN;
      rtemp[i] = comoving_distance(ztemp[i]);
  }
  red_acc = gsl_interp_accel_alloc(), dist_acc = gsl_interp_accel_alloc();
  red_spline = gsl_spline_alloc(gsl_interp_cspline, nbins), dist_spline = gsl_spline_alloc(gsl_interp_cspline, nbins);
  gsl_spline_init(red_spline, rtemp, ztemp, nbins), gsl_spline_init(dist_spline, ztemp, rtemp, nbins);
  free(rtemp), free(ztemp);

  return;
}

// Destroy the grid(s)
// ===================
void destroy_grids(void) {
  fftw_destroy_plan(plan);
  free(ddg);
  if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) {
  	fftw_destroy_plan(plan_mom);
  	free(ddg_mom);
  }
  if (Survey) {
    free(ddg_2);
    free(F0);
    if ((Momentum != 0) && (Momentum != 1)) {
      free(ddg_mom_2);
      free(F0_mom);
    }
  }
  if (DoInterlacing) {
    fftw_destroy_plan(plan_interlace);
    free(ddg_interlace);
    if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) {
  	  fftw_destroy_plan(plan_mom_interlace);
  	  free(ddg_mom_interlace);
    }
    if (Survey) {
      fftw_destroy_plan(plan_interlace);
      free(ddg_interlace_2);
      if ((Momentum != 0) && (Momentum != 1)) {
        free(ddg_mom_interlace_2);
      }
    }
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
