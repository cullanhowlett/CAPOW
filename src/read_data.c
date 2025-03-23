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

/* ====================================================================================*/
/* This file contains routines to read in data in ASCII or GADGET format.              */
/* Contains serial and parallel versions and also routines to assign data to the grid. */
/* ====================================================================================*/

#include "vars.h"
#include "proto.h"

// Read in an ASCII file containing the simulation box data. 
// Every processor reads in the file but only stores the relevant parts.
// =====================================================================
double read_periodic_serial_ascii(char *inputfile) {

  FILE * fp;
  int bufsize = 2000;
  char buf[bufsize];
  double NREAD = 0, NOUT = 0, NGRID = 0, NPV = 0;
  double XMIN_LOCAL = Local_x_start*dx+XMIN;
  double XMAX_LOCAL = (Local_x_start+Local_nx)*dx+XMIN;

  if (ThisTask == 0) {
    printf("\nReading: %s\n", inputfile); 
    fflush(stdout);
  }

  if (Momentum) {
    nsq = 0.0; vr_ave = 0.0; vrsq_ave = 0.0; vr3_ave = 0.0; vr4_ave = 0.0;
  }

  // Open the file
  if((fp=fopen(inputfile,"r"))==NULL) { printf("Task %d cannot open input file\n", ThisTask);  FatalError("read_data", 96); }

  // Loop over each line and assign the particles to the grid
  // We skip over any lines starting with `#' as we assume they are headers
  // The exact format will likely need modifying
  while(fgets(buf,bufsize,fp)) {
    if(strncmp(buf,"#",1)==0) continue;

    // These lines will probably need changing to match the input data format.
    double w = 1.0;
    double tx,ty,tz,tvx,tvy,tvz,tvr;
    //if(sscanf(buf,"%lf %lf %lf %lf\n",&tx,&ty,&tz,&tz_rsd)!=4) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
    //tz = tz_rsd;
    if(sscanf(buf,"%lf %lf %lf %lf %lf %lf\n",&tx,&ty,&tz,&tvx,&tvy,&tvz)!=6) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
    if (LOS == 0) {
      tvr = tvz;
    } else if (LOS == 1) {
      tvr = tvx;
      tx += tvr * (1.0 + Redshift) / (100.0*sqrt(Omega_m*pow((1.0 + Redshift),3.0) + 1.0 - Omega_m));
    } else if (LOS == 2) {
      tvr = tvy;
      ty += tvr * (1.0 + Redshift) / (100.0*sqrt(Omega_m*pow((1.0 + Redshift),3.0) + 1.0 - Omega_m));
    } else {
      tvr = tvz;
      tz += tvr * (1.0 + Redshift) / (100.0*sqrt(Omega_m*pow((1.0 + Redshift),3.0) + 1.0 - Omega_m));
    }

    NREAD += 1;

    // If we are adding RSD we'll just highlight how many fall outside the boundaries, but won't error out.
    if ((tx < XMIN) || (tx > XMAX) || (ty < YMIN) || (ty > YMAX) || (tz < ZMIN) || (tz > ZMAX)) {
      if (LOS == 0) {
        printf("Task %d has particle out of grid bounds: x=%lf, y=%lf, z=%lf\n", ThisTask, tx, ty, tz);
        if (LOS == 0) FatalError("read_data", 110);
      } else {
        NOUT += 1;
      }
    }

    // P00, P01, P11, P02, P12, P22, P03, P13, P04
    if (Momentum == 0) {
      NGRID += add_to_grid(tx, ty, tz, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace); 
    } else if (Momentum == 1) {
      NGRID += add_to_grid(tx, ty, tz, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace); 
      NPV += add_to_grid(tx, ty, tz, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);
    } else if (Momentum == 2) {
      NGRID += add_to_grid(tx, ty, tz, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
    } else if (Momentum == 3) {
      NGRID += add_to_grid(tx, ty, tz, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
      NPV += add_to_grid(tx, ty, tz, w, tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);
    } else if (Momentum == 4) {
      NGRID += add_to_grid(tx, ty, tz, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
      NPV += add_to_grid(tx, ty, tz, w, tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);
    } else if (Momentum == 5) {
      NGRID += add_to_grid(tx, ty, tz, w, tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
    } else if (Momentum == 6) {
      NGRID += add_to_grid(tx, ty, tz, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
      NPV += add_to_grid(tx, ty, tz, w, tvr*tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr*tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);
    } else if (Momentum == 7) {
      NGRID += add_to_grid(tx, ty, tz, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
      NPV += add_to_grid(tx, ty, tz, w, tvr*tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr*tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);
    } else {
      NGRID += add_to_grid(tx, ty, tz, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
      NPV += add_to_grid(tx, ty, tz, w, tvr*tvr*tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
      if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, tvr*tvr*tvr*tvr, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);
    }

    if (Momentum) {
      nsq += w;
      vr_ave += w*tvr;
      vrsq_ave += w*tvr*tvr;
      vr3_ave += w*tvr*tvr*tvr;
      vr4_ave += w*tvr*tvr*tvr*tvr;
    }

  }

  if (LOS == 0) {
    printf("Task %d read %lf particles, gridded %lf\n", ThisTask, NREAD, NGRID);
    fflush(stdout);
  } else {
    printf("Task %d read %lf particles, gridded %lf (%lf out of bounds due to RSD)\n", ThisTask, NREAD, NGRID, NOUT);
    fflush(stdout);
  }

  // Copy across the extra slices from the task on the left and add it to the leftmost slices
  // of the task on the right. Skip over tasks without any slices.
  if (InterpOrder > 0) {
    double * temp_ddg = (double *)calloc(InterpOrder*alloc_slice,sizeof(double));
    ierr = MPI_Sendrecv(&(ddg[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                        &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
    for (int i=0;i<InterpOrder*alloc_slice;i++) ddg[i] += temp_ddg[i];
    if (DoInterlacing) {
      ierr = MPI_Sendrecv(&(ddg_interlace[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                          &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
      for (int i=0;i<InterpOrder*alloc_slice;i++) ddg_interlace[i] += temp_ddg[i];
    }
    free(temp_ddg);
    if ((Momentum != 0) && (Momentum != 2) && (Momentum != 5)) {
      double * temp_ddg = (double *)calloc(InterpOrder*alloc_slice,sizeof(double));
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

  return NGRID;

}

// Read in an ASCII file containing the survey data.
// Every processor reads in the file but only stores the relevant parts.
// =====================================================================
double read_survey_serial_ascii(char *inputfile, struct survey_data * inputdata, int randoms) {

  FILE * fp;
  int bufsize = 2000;
  int largebufsize = 8000000;
  char largebuf[largebufsize];
  unsigned long long NREAD = 0, NKEEP = 0;

  double XLOW = 1.0e-30, XHI = -1.0e30; 
  double YLOW = 1.0e-30, YHI = -1.0e30;
  double ZLOW = 1.0e-30, ZHI = -1.0e30;

  double XMIN_Interp = XMIN + (InterpOrder-1)*dx, XMAX_Interp = XMAX - (InterpOrder-1)*dx;
  double YMIN_Interp = YMIN + (InterpOrder-1)*dy, YMAX_Interp = YMAX - (InterpOrder-1)*dy;
  double ZMIN_Interp = ZMIN + (InterpOrder-1)*dz, ZMAX_Interp = ZMAX - (InterpOrder-1)*dz;
  if (DoInterlacing) {
    XMAX_Interp -= 0.5*dx;
    YMAX_Interp -= 0.5*dy;
    ZMAX_Interp -= 0.5*dz;
  }

  if (ThisTask == 0) {
    printf("Reading: %s\n", inputfile); 
    fflush(stdout);
  }

  // Open the file
  if((fp=fopen(inputfile,"r"))==NULL) { printf("Task %d cannot open input file\n", ThisTask);  FatalError("read_data", 113); }

  // Skip over lines starting with #
  while(fgets(largebuf,largebufsize,fp)) {
    if(strncmp(largebuf,"#",1)!=0) break;
  }

  // Move back to the start of the line
  fseek(fp, -strlen(largebuf), SEEK_CUR);

  // Loop over each line and store the data coordinates
  // The exact format may need modifying
  int nleft=0;
  char leftover[bufsize];

  int nbuf = fread(largebuf+nleft, 1, largebufsize-nleft-1, fp);
  largebuf[largebufsize-1] = '\0';
  //printf("%d\n", nbuf + nleft);
  do {

    //printf("%s\n", largebuf);

    int tid, tcen;
    double tx, ty, tz, tred, tnbar, tmass, tw, tdist;
    double tlogdist, tlogdist_true, tlogdist_err;

    char * buf = strtok(largebuf, "\n");
    unsigned long largelen = strlen(buf)+1;
    while(largelen <= nbuf+nleft) {
      if (randoms) {
        //if(sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf\n",&tx,&ty,&tz,&tlogdist,&tlogdist_true,&tlogdist_err,&tnbar)!=7) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
        if(sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&tx,&ty,&tz,&tlogdist_true,&tlogdist_err,&tlogdist,&tmass,&tmass,&tmass,&tnbar)!=10) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
        tw = 1.0;
      } else {
        //if(sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf\n",&tx,&ty,&tz,&tlogdist,&tlogdist_true,&tlogdist_err,&tnbar)!=7) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
        if(sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&tx,&ty,&tz,&tlogdist_true,&tlogdist_err,&tlogdist,&tmass,&tmass,&tmass,&tnbar)!=10) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
        tw = 1.0;
      }
      NREAD++;
      
      if (ThisTask == 0) {
        if ((NREAD % 1000000) == 0) {
          printf("Read %llu objects, kept %llu\n", NREAD, NKEEP);
          fflush(stdout);
        }
      }

      if (Coord_Type == 0) {
        tred = gsl_spline_eval(red_spline, sqrt(tx*tx + ty*ty + tz*tz), red_acc);
        if ((tred < REDMININ) || (tred > REDMAXIN)) {
          buf = strtok(NULL, "\n");
          if (buf == NULL) break;
          largelen += strlen(buf)+1;
          continue;
        }
      } else {
        tred = tz;
        if ((tred < REDMININ) || (tred > REDMAXIN)) {
          buf = strtok(NULL, "\n");
          if (buf == NULL) break;
          largelen += strlen(buf)+1;
          continue;
        }
        //double dist = comoving_distance(tred);
        tdist = gsl_spline_eval(dist_spline, tred, dist_acc);
        double ra = tx, dec = ty;
        if (Coord_Type == 1) {
          ra *= (M_PI/180.0);
          dec *= (M_PI/180.0);
        }
        tx = tdist*cos(dec)*cos(ra);
        ty = tdist*cos(dec)*sin(ra);
        tz = tdist*sin(dec);
        if (Momentum && randoms) tlogdist_err = 0.1;
      }

      // Check the data can fit in the grid including whatever interpolation order we are using
      if ((tx < XMIN_Interp) || (tx >= XMAX_Interp) || (ty < YMIN_Interp) || (ty >= YMAX_Interp) || (tz < ZMIN_Interp) || (tz >= ZMAX_Interp)) {
        printf("Task %d has object out of grid bounds for chosen InterpOrder: x=%lf, y=%lf, z=%lf\n", ThisTask, tx, ty, tz);
        FatalError("read_data", 155);
      }

      // Add it to the data structure
      inputdata[NKEEP].coord[0] = tx;
      inputdata[NKEEP].coord[1] = ty;
      inputdata[NKEEP].coord[2] = tz;
      inputdata[NKEEP].redshift = tred;
      inputdata[NKEEP].weight = tw;
      if (NBAR_Column > 0) inputdata[NKEEP].nbar = tnbar;
      if (Momentum) {
        double q = 3.0/2.0*Omega_m - 1.0;
        double zmod = tred*(1.0 + 0.5*(1.0 - q)*tred - 1.0/6.0*(2.0 - q - 3.0*q*q)*tred*tred);
        inputdata[NKEEP].pv = LightSpeed*zmod*log(10.0)*tlogdist/(1.0+zmod);
        if (randoms) {
          inputdata[NKEEP].pverr = LightSpeed*zmod*log(10.0)*tlogdist_err/(1.0+zmod);
        } else {
          inputdata[NKEEP].pverr = LightSpeed*zmod*log(10.0)*tlogdist_err/(1.0+zmod);
        }
      }
      NKEEP++;

      if (tred < REDMIN) REDMIN = tred;
      if (tred > REDMAX) REDMAX = tred;

      if (tx < XLOW) XLOW = tx;
      if (ty < YLOW) YLOW = ty;
      if (tz < ZLOW) ZLOW = tz;
      if (tx > XHI) XHI = tx;
      if (ty > YHI) YHI = ty;
      if (tz > ZHI) ZHI = tz;

      buf = strtok(NULL, "\n");
      if (buf == NULL) break;
      largelen += strlen(buf)+1;
    }

    if (feof(fp)) break;

    if (buf == NULL) {
      nleft = 0;
    } else {
      nleft = strlen(buf);
      memcpy(leftover, buf, nleft);
      leftover[nleft] = '\0';
      memcpy(largebuf, leftover, nleft);
    }
    nbuf = fread(largebuf+nleft, 1, largebufsize-nleft-1, fp);
    largebuf[largebufsize-1] = '\0';
    //printf("%d\n", nbuf + nleft);

  } while (!ferror(fp));

  if (ThisTask == 0) {
    printf("Have %llu objects\n", NKEEP);
    printf("%12.6lf <   X  < %12.6lf\n", XLOW, XHI);
    printf("%12.6lf <   Y  < %12.6lf\n", YLOW, YHI);
    printf("%12.6lf <   Z  < %12.6lf\n", ZLOW, ZHI);
  }

  fclose(fp);

  return (double)NKEEP;
}

// Calculates the comoving distance from the redshift
// ==================================================
double comoving_distance(double red) {
  double result, error;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  F.function = &f;
  gsl_integration_qags(&F, 0.0, red, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return LightSpeed*result/100.0;
}

// Integrand for the comoving distance
// ===================================
double f(double z, void *p) {
  double ff = 1.0/sqrt(Omega_m*(1.0+z)*(1.0+z)*(1.0+z)+(1.0-Omega_m));
  return ff;
}

// Compute the number density of the data. Has a flag to tell it 
// whether each processor contains all the data, or a unique portion.
// ==================================================================
void compute_nbar(int parallel, unsigned long long NDATA, unsigned long long NRAND) {

  if (ThisTask == 0) printf("Calculating Number Density...\n");

  int nbins = (int)(ceil((REDMAXIN-REDMININ)/0.002));
  double * znbar = (double *)calloc(nbins+2, sizeof(double));
  double * nbar = (double *)calloc(nbins+2, sizeof(double));

  if (parallel) { 
    double REDMIN_glob, REDMAX_glob;
    MPI_Allreduce(&REDMIN, &REDMIN_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&REDMAX, &REDMAX_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    REDMIN = REDMIN_glob;
    REDMAX = REDMAX_glob;
  }
  double redbinwidth = (REDMAX-REDMIN)/(float)nbins;

  // Compute the number density in redshift bins
  for (unsigned long long i=0; i<NDATA; i++) {
    int nbin = (int)(floor((randoms[i].redshift-REDMIN)/redbinwidth))+1;
    if (nbin == nbins+1) nbin--;
    nbar[nbin]++;
  }

  for (int i=0;i<nbins;i++) {
    double outervolume = gsl_spline_eval(dist_spline, (i+1)*redbinwidth+REDMIN, dist_acc);
    double innervolume = gsl_spline_eval(dist_spline, i*redbinwidth+REDMIN, dist_acc);
    double volume = SkyArea/(180.0*180.0)*M_PI*M_PI*(outervolume*outervolume*outervolume-innervolume*innervolume*innervolume)/3.0;
    znbar[i+1] = (i+0.5)*redbinwidth+REDMIN;
    nbar[i+1] /= volume;
    printf("%12.6lf, %12.6lf\n", znbar[i+1], nbar[i+1]);
  }

  if (parallel) { 
    double * nbar_glob = (double *)calloc(nbins+2, sizeof(double));
    MPI_Allreduce(nbar, nbar_glob, nbins+2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   for (int i=0; i<nbins+2; i++) nbar[i] = nbar_glob[i];
  }

  // Ensure the number density goes to zero at the endpoints.
  znbar[0] = REDMIN-0.5*redbinwidth;
  znbar[nbins+1] = REDMAX+0.5*redbinwidth;

  // Spline interpolate the number density
  gsl_interp_accel * nbar_acc = gsl_interp_accel_alloc();
  gsl_spline * nbar_spline = gsl_spline_alloc(gsl_interp_cspline, nbins+2);
  gsl_spline_init(nbar_spline, znbar, nbar, nbins+2);

  // Assign the number density to each data and random point
  for (unsigned long long i=0; i<NDATA; i++) data[i].nbar = gsl_spline_eval(nbar_spline, data[i].redshift, nbar_acc);
  for (unsigned long long i=0; i<NRAND; i++) randoms[i].nbar = gsl_spline_eval(nbar_spline, randoms[i].redshift, nbar_acc);
  //if (Momentum != 1) {
  //  for (unsigned long long i=0; i<NRAND; i++) randoms[i].nbar = gsl_spline_eval(nbar_spline, randoms[i].redshift, nbar_acc);
  //}

  free(znbar);
  free(nbar);
  gsl_spline_free(nbar_spline);
  gsl_interp_accel_free(nbar_acc);

  return;
}

// Compute the FKP weights for each object
// =======================================
void compute_fkp(unsigned long long NOBJ, struct survey_data * inputdata) {
  for (unsigned long long i=0; i<NOBJ; i++) {
    if (Momentum == 0) {
      inputdata[i].weight *= 1.0/(1.0+inputdata[i].nbar*FKP_Pk);
    } else if (Momentum == 1) {
      inputdata[i].weight *= 1.0/((inputdata[i].pverr*inputdata[i].pverr+300.0*300.0)+inputdata[i].nbar*FKP_Pk_mom);
    } else {
      double temp_weight = 1.0/(1.0+inputdata[i].nbar*FKP_Pk);
      double temp_weight_pv = 1.0/((inputdata[i].pverr*inputdata[i].pverr+300.0*300.0)+inputdata[i].nbar*FKP_Pk_mom);
      double temp_weight_cross = 1.0/sqrt(1.0 + temp_weight*temp_weight_pv*inputdata[i].nbar*inputdata[i].nbar*FKP_Pk_cross*FKP_Pk_cross);
      inputdata[i].weight_pv = inputdata[i].weight*temp_weight_pv*temp_weight_cross;
      inputdata[i].weight *= temp_weight*temp_weight_cross;
      //printf("%g, %g, %g\n", temp_weight_pv, inputdata[i].weight, inputdata[i].weight_pv);
    }
  }
  return;
}

// Assign survey data to a grid. The third argument is the factor to multiply
// each assignment by which allows the function to be used for both data and randoms.
// ==================================================================================
double assign_survey_data(unsigned long long NOBJ, struct survey_data * inputdata, double prefactor, int data) {

  double NGRID = 0;
  double XMIN_LOCAL = Local_x_start*dx+XMIN;
  double XMAX_LOCAL = (Local_x_start+Local_nx)*dx+XMIN;

  for (unsigned long long i=0; i<NOBJ; i++) {
    if (Momentum == 1) {
      NGRID += add_to_grid(inputdata[i].coord[0], inputdata[i].coord[1], inputdata[i].coord[2], inputdata[i].weight, prefactor*inputdata[i].pv, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg)/(inputdata[i].pv*inputdata[i].weight);
      if (DoInterlacing) add_to_grid(inputdata[i].coord[0]+dx/2.0, inputdata[i].coord[1]+dy/2.0, inputdata[i].coord[2]+dz/2.0, inputdata[i].weight, prefactor*inputdata[i].pv, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
    } else {
      NGRID += add_to_grid(inputdata[i].coord[0], inputdata[i].coord[1], inputdata[i].coord[2], prefactor*inputdata[i].weight, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg)/inputdata[i].weight;
      if (DoInterlacing) add_to_grid(inputdata[i].coord[0]+dx/2.0, inputdata[i].coord[1]+dy/2.0, inputdata[i].coord[2]+dz/2.0, prefactor*inputdata[i].weight, 1.0, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
      if (Momentum && data) {
        add_to_grid(inputdata[i].coord[0], inputdata[i].coord[1], inputdata[i].coord[2], inputdata[i].weight_pv, prefactor*inputdata[i].pv, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom);
        if (DoInterlacing) add_to_grid(inputdata[i].coord[0]+dx/2.0, inputdata[i].coord[1]+dy/2.0, inputdata[i].coord[2]+dz/2.0, inputdata[i].weight_pv, prefactor*inputdata[i].pv, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_mom_interlace);  
      }
    }
  }
  printf("Task %d has %llu objects, gridded %lf\n", ThisTask, NOBJ, NGRID/prefactor);
  fflush(stdout);

  return NGRID;
}

// Read in multiple ASCII files containing the data. 
// Every processor reads in a unique file and the processors sort out where
// it needs to go to be allocated to a grid before transferring approriately.
// ==========================================================================
double read_periodic_parallel_ascii(void) {

	if (ThisTask == 0) printf("\nERROR: read_periodic_parallel_ascii not implemented yet.\n");
  FatalError((char *)"read_param.c", 378);

  return 0.0;

}


// Add a particle to the grid using the given order of interpolation. Wraps the particles for a periodic box
// Written in such a way we only have to create and copy buffer slices of the grid on the right hand edge.
// For survey data there is no wrapping, but we have already checked that there is extra grid slices on the outside to assign to.
// ==============================================================================================================================
double add_to_grid(double x, double y, double z, double w, double vr, double xmin, double xmax, int nx, double * density) {

  double wvr = w * vr;

  // Nearest Neighbour interpolation
  if (InterpOrder == 1) {

    if (ThisTask == NTask-1) {
      if ((x < xmin) || (x >= xmax)) return 0.0;
    } else { 
      if ((x < xmin) || (x > xmax)) return 0.0;
    }
    if ((y < YMIN) || (y >= YMAX)) return 0.0;
    if ((z < ZMIN) || (z >= ZMAX)) return 0.0;

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
    density[ind] += wvr;

  // Cloud-in-Cell interpolation
  } else if (InterpOrder == 2) {

    if (ThisTask == NTask-1) {
      if ((x < xmin) || (x >= xmax)) return 0.0;
    } else { 
      if ((x < xmin) || (x > xmax)) return 0.0;
    }
    if ((y < YMIN) || (y >= YMAX)) return 0.0;
    if ((z < ZMIN) || (z >= ZMAX)) return 0.0;

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
    if (Periodic) {
      if(ixneigh == nx) ixneigh = nx-1;
      if (LOS == 0) {
        if(iyneigh == NY) iyneigh = 0;
        if(izneigh == NZ) izneigh = 0;
      }
    } else if (Survey) {
      if(ixneigh == nx) {
        ixneigh = nx-1;
        ix = nx-2;
      }
      if(iyneigh == NY) {
        iyneigh = NY-1;
        iy = NY-2;
      }
      if(izneigh == NZ) {
        izneigh = NZ-1;
        iz = NZ-2;
      }
    }

    density[iz+2*(NZ/2+1)*(iy+NY*ix)]           += wvr*itx*ity*itz;
    density[iz+2*(NZ/2+1)*(iy+NY*ixneigh)]      += wvr*idx*ity*itz;
    if (iyneigh < NY) {
      density[iz+2*(NZ/2+1)*(iyneigh+NY*ix)]      += wvr*itx*idy*itz;
      density[iz+2*(NZ/2+1)*(iyneigh+NY*ixneigh)] += wvr*idx*idy*itz;
      if (izneigh < NZ) {
        density[izneigh+2*(NZ/2+1)*(iyneigh+NY*ix)]      += wvr*itx*idy*idz;
        density[izneigh+2*(NZ/2+1)*(iyneigh+NY*ixneigh)] += wvr*idx*idy*idz;
      }
    }
    if (izneigh < NZ) {
      density[izneigh+2*(NZ/2+1)*(iy+NY*ix)]      += wvr*itx*ity*idz;
      density[izneigh+2*(NZ/2+1)*(iy+NY*ixneigh)] += wvr*idx*ity*idz;
    }

  // Triangular-Shaped Cloud interpolation
  } else if (InterpOrder == 3) {

    if ((Periodic) && (LOS == 0)) {
      if (x-XMIN < dx) x += XMAX-XMIN;
      if (ThisTask == NTask-1) {
        if ((x < xmin+dx) || (x >= xmax+dx)) return 0.0;
      } else { 
        if ((x < xmin+dx) || (x > xmax+dx)) return 0.0;
      }
    } else {
      if (ThisTask == 0) {
        if ((x < xmin) || (x >= xmax+dx)) return 0.0;
      } else if (ThisTask == NTask-1) {
        if ((x < xmin+dx) || (x >= xmax)) return 0.0;
      } else { 
        if ((x < xmin+dx) || (x > xmax+dx)) return 0.0;
      }
    }
    if ((y < YMIN) || (y >= YMAX)) return 0.0;
    if ((z < ZMIN) || (z >= ZMAX)) return 0.0;

    double scalex = (x-xmin)/dx;
    double scaley = (y-YMIN)/dy;
    double scalez = (z-ZMIN)/dz;

    int ix = (int)(scalex+0.5);
    int iy = (int)(scaley+0.5);
    int iz = (int)(scalez+0.5);
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

    if (Periodic) {
      if(ixneighhi == nx) ixneighhi = nx-1;
      if (LOS == 0) {
        if (iyneighlow < 0) iyneighlow = NY-1;
        if (izneighlow < 0) izneighlow = NZ-1;
        if(iyneighhi >= NY) {
          iyneighhi -= NY;
          if(iy >= NY) iy = 0;
          if(iyneighhi > 1) iyneighhi = 1;
        }
        if(izneighhi >= NZ) {
          izneighhi -= NZ;
          if(iz >= NZ) iz = 0;
          if(izneighhi > 1) izneighhi = 1;
        }
      }
    } else if (Survey) {
      if(ixneighhi == nx) {
        ixneighhi = nx-1;
        ix = nx-2;
        ixneighlow = nx-3;
      }
      if(iyneighhi == NY) {
        iyneighhi = NY-1;
        iy = NY-2;
        iyneighlow = NY-3;
      }
      if(izneighhi == NZ) {
        izneighhi = NZ-1;
        iz = NZ-2;
        izneighlow = NZ-3;
      }
    }

    if ((iy < NY) && (iz < NZ)) {
      density[iz+2*(NZ/2+1)*(iy+NY*ix)]                             += wvr*itx*ity*itz;
      if (ixneighlow > 0) density[iz+2*(NZ/2+1)*(iy+NY*ixneighlow)] += wvr*isx*ity*itz;
      if (ixneighhi < nx) density[iz+2*(NZ/2+1)*(iy+NY*ixneighhi)]  += wvr*idx*ity*itz;
    }
    if (iyneighlow > 0) {
      if (iz < NZ) {
        density[iz+2*(NZ/2+1)*(iyneighlow+NY*ix)]                             += wvr*itx*isy*itz;
        if (ixneighlow > 0) density[iz+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)] += wvr*isx*isy*itz;
        if (ixneighhi < nx) density[iz+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]  += wvr*idx*isy*itz;
      }
      if (izneighlow > 0) {
            density[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ix)]                             += wvr*itx*isy*isz;
            if (ixneighlow > 0) density[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)] += wvr*isx*isy*isz;
            if (ixneighhi < nx) density[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]  += wvr*idx*isy*isz;
      }
      if (izneighhi < NZ) {
            density[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ix)]                             += wvr*itx*isy*idz;
            if (ixneighlow > 0) density[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)] += wvr*isx*isy*idz;
            if (ixneighhi < nx) density[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]  += wvr*idx*isy*idz;
      }
    }
    if (iyneighhi < NY) {
      if (iz < NZ) {
        density[iz+2*(NZ/2+1)*(iyneighhi+NY*ix)]                             += wvr*itx*idy*itz;
        if (ixneighlow > 0) density[iz+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)] += wvr*isx*idy*itz;
        if (ixneighhi < nx) density[iz+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]  += wvr*idx*idy*itz;
      }
      if (izneighlow > 0) {
        density[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ix)]                             += wvr*itx*idy*isz;
        if (ixneighlow > 0) density[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)] += wvr*isx*idy*isz;
        if (ixneighhi < nx) density[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]  += wvr*idx*idy*isz;
      }
      if (izneighhi < NZ) {
        density[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ix)]                             += wvr*itx*idy*idz;
        if (ixneighlow > 0) density[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)] += wvr*isx*idy*idz;
        if (ixneighhi < nx) density[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]  += wvr*idx*idy*idz;
      }
    }
    if (izneighlow > 0) {
      if (iy < NY) {
        density[izneighlow+2*(NZ/2+1)*(iy+NY*ix)]                             += wvr*itx*ity*isz;
        if (ixneighlow > 0) density[izneighlow+2*(NZ/2+1)*(iy+NY*ixneighlow)] += wvr*isx*ity*isz;
        if (ixneighhi < nx) density[izneighlow+2*(NZ/2+1)*(iy+NY*ixneighhi)]  += wvr*idx*ity*isz;
      }
    }
    if (izneighhi < NZ) {
      if (iy < NY) {
        density[izneighhi+2*(NZ/2+1)*(iy+NY*ix)]                             += wvr*itx*ity*idz;
        if (ixneighlow > 0) density[izneighhi+2*(NZ/2+1)*(iy+NY*ixneighlow)] += wvr*isx*ity*idz;
        if (ixneighhi < nx) density[izneighhi+2*(NZ/2+1)*(iy+NY*ixneighhi)]  += wvr*idx*ity*idz;
      }
    }
  }
  return w;
}
