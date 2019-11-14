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
double read_periodic_serial_ascii(void) {

  FILE * fp;
  int bufsize = 2000;
  char buf[bufsize], inputfile[bufsize];
  double NREAD = 0, NGRID = 0;
  double XMIN_LOCAL = Local_x_start*dx+XMIN;
  double XMAX_LOCAL = (Local_x_start+Local_nx)*dx+XMIN;
  sprintf(inputfile, "%s/%s", InputDir, FileBase);

  if (ThisTask == 0) {
    printf("\nReading: %s\n", inputfile); 
    fflush(stdout);
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
    double tx,ty,tz,tz_rsd;
    if(sscanf(buf,"%lf %lf %lf %lf\n",&tx,&ty,&tz,&tz_rsd)!=4) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }
    tz = tz_rsd;

    NREAD += w;
    if ((tx < XMIN) || (tx >= XMAX) || (ty < YMIN) || (ty >= YMAX) || (tz < ZMIN) || (tz >= ZMAX)) {
      printf("Task %d has particle out of grid bounds: x=%lf, y=%lf, z=%lf\n", ThisTask, tx, ty, tz);
      FatalError("read_data", 110);
    }

    NGRID += add_to_grid(tx, ty, tz, w, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
    if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);

  }
  printf("Task %d read %lf particles, gridded %lf\n", ThisTask, NREAD, NGRID);
  fflush(stdout);

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
  }

  return NGRID;

}

// Read in an ASCII file containing the survey data.
// Every processor reads in the file but only stores the relevant parts.
// =====================================================================
double read_survey_serial_ascii(char *inputfile, struct survey_data * inputdata) {

  FILE * fp;
  int bufsize = 2000;
  char buf[bufsize];
  unsigned long long NREAD = 0;

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

  // Loop over each line and store the data coordinates
  // We skip over any lines starting with `#' as we assume they are headers
  // The exact format may need modifying
  int strlen = x_Column;
  if (y_Column > strlen) strlen = y_Column;
  if (z_Column > strlen) strlen = z_Column;
  if (NBAR_Column > strlen) strlen = NBAR_Column;
  if (FKP_Column > strlen) strlen = FKP_Column;
  while(fgets(buf,bufsize,fp)) {
    if(strncmp(buf,"#",1)==0) continue;

    char *bufp = buf;
    int offset;
    double val;
    double tx, ty, tz, tred, tnbar, tw;
    for (int i=1; i<=strlen; i++) {
      sscanf(bufp, "%lf%n", &val, &offset);
      if (i == x_Column) tx = val;
      if (i == y_Column) ty = val;
      if (i == z_Column) {
        tz = val;
        if (Coord_Type != 0) tred = tz;
      }
      if (i == NBAR_Column) tnbar = val;
      if (i == FKP_Column) tw = val;
      bufp += offset;
    }
    //double tx, ty, tz, tz_rsd, tred, tnbar, tw;
    //if(sscanf(buf,"%lf %lf %lf %lf\n",&tx,&ty,&tz,&tred)!=4) { printf("Task %d has error reading file: %s\n", ThisTask, buf);  FatalError("read_data", 102); }

    // Compute the redshift if we only had cartesian coordinates, otherwise convert to cartesian coordinates
    if (Coord_Type == 0) {
      tred = gsl_spline_eval(red_spline, sqrt(tx*tx + ty*ty + tz*tz), red_acc);
    } else {
      //double dist = comoving_distance(tred);
      double dist = gsl_spline_eval(dist_spline, tred, dist_acc);
      double ra = tx, dec = ty;
      if (Coord_Type == 1) {
        ra *= (M_PI/180.0);
        dec *= (M_PI/180.0);
      }
      tx = dist*cos(dec)*cos(ra);
      ty = dist*cos(dec)*sin(ra);
      tz = dist*sin(dec);
    }

    // Check the data can fit in the grid including whatever interpolation order we are using
    if ((tx < XMIN_Interp) || (tx >= XMAX_Interp) || (ty < YMIN_Interp) || (ty >= YMAX_Interp) || (tz < ZMIN_Interp) || (tz >= ZMAX_Interp)) {
      printf("Task %d has object out of grid bounds for chosen InterpOrder: x=%lf, y=%lf, z=%lf\n", ThisTask, tx, ty, tz);
      FatalError("read_data", 155);
    }

    // Add it to the data structure
    inputdata[NREAD].coord[0] = tx;
    inputdata[NREAD].coord[1] = ty;
    inputdata[NREAD].coord[2] = tz;
    inputdata[NREAD].redshift = tred;
    if (NBAR_Column > 0) inputdata[NREAD].nbar = tnbar;
    if (FKP_Column > 0) inputdata[NREAD].weight = tw;
    NREAD ++;

    if (tred < REDMIN) REDMIN = tred;
    if (tred > REDMAX) REDMAX = tred;

    if (tx < XLOW) XLOW = tx;
    if (ty < YLOW) YLOW = ty;
    if (tz < ZLOW) ZLOW = tz;
    if (tx > XHI) XHI = tx;
    if (ty > YHI) YHI = ty;
    if (tz > ZHI) ZHI = tz;

    if (ThisTask == 0) {
      if ((NREAD % 1000000) == 0) {
        printf("Read %llu objects\n", NREAD);
        fflush(stdout);
      }
    }

  }

  if (ThisTask == 0) {
    printf("Read %llu objects\n", NREAD);
    printf("%12.6lf <   X  < %12.6lf\n", XLOW, XHI);
    printf("%12.6lf <   Y  < %12.6lf\n", YLOW, YHI);
    printf("%12.6lf <   Z  < %12.6lf\n", ZLOW, ZHI);
  }

  return (double)NREAD;
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

  printf("Calculating Number Density...\n");

  int nbins = 200;
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
    int nbin = (int)(floor((data[i].redshift-REDMIN)/redbinwidth))+1;
    if (nbin == nbins+1) nbin--;
    nbar[nbin]++;
  }

  for (int i=0;i<nbins;i++) {
    double outervolume = gsl_spline_eval(dist_spline, (i+1)*redbinwidth+REDMIN, dist_acc);
    double innervolume = gsl_spline_eval(dist_spline, i*redbinwidth+REDMIN, dist_acc);
    double volume = SkyArea/(180.0*180.0)*M_PI*M_PI*(outervolume*outervolume*outervolume-innervolume*innervolume*innervolume)/3.0;
    znbar[i+1] = (i+0.5)*redbinwidth+REDMIN;
    nbar[i+1] /= volume;
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
  for (unsigned long long i=0; i<NDATA; i++) {
    data[i].nbar = gsl_spline_eval(nbar_spline, data[i].redshift, nbar_acc);
    //if (fabs(data[i].nbar) > 1.0e-3) printf("%lf, %lf\n", data[i].redshift, data[i].nbar);
  }
  for (unsigned long long i=0; i<NRAND; i++) randoms[i].nbar = gsl_spline_eval(nbar_spline, randoms[i].redshift, nbar_acc);

  free(znbar);
  free(nbar);
  gsl_spline_free(nbar_spline);
  gsl_interp_accel_free(nbar_acc);

  return;
}

// Compute the FKP weights for each object
// =======================================
void compute_fkp(unsigned long long NOBJ, struct survey_data * inputdata) {
  for (unsigned long long i=0; i<NOBJ; i++) inputdata[i].weight = 1.0/(1.0+inputdata[i].nbar*FKP_Pk);
  return;
}

// Assign survey data to a grid. The third argument is the factor to multiply
// each assignment by which allows the function to be used for both data and randoms.
// ==================================================================================
double assign_survey_data(unsigned long long NOBJ, struct survey_data * inputdata, double prefactor) {

  double NGRID = 0;
  double XMIN_LOCAL = Local_x_start*dx+XMIN;
  double XMAX_LOCAL = (Local_x_start+Local_nx)*dx+XMIN;

  for (unsigned long long i=0; i<NOBJ; i++) {
    NGRID += add_to_grid(inputdata[i].coord[0], inputdata[i].coord[1], inputdata[i].coord[2], prefactor*inputdata[i].weight, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg);
    if (DoInterlacing) add_to_grid(inputdata[i].coord[0]+dx/2.0, inputdata[i].coord[1]+dy/2.0, inputdata[i].coord[2]+dz/2.0, prefactor*inputdata[i].weight, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_interlace);
  }
  printf("Task %d has %llu objects, gridded %lf\n", ThisTask, NOBJ, NGRID/prefactor);
  fflush(stdout);

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
  }

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
double add_to_grid(double x, double y, double z, double w, double xmin, double xmax, int nx, double * density) {

  // Nearest Neighbour interpolation
  if (InterpOrder == 1) {

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
    density[ind] += w;

  // Cloud-in-Cell interpolation
  } else if (InterpOrder == 2) {

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
    if (Periodic) {
      if(ixneigh >= nx) ixneigh = nx;
      if(iyneigh >= NY) {
        iyneigh = 0;
        if (iy >= NY) iy = NY-1;
      }
      if(izneigh >= NZ) {
        izneigh = 0;
        if (iz >= NZ) iz = NZ-1;
      }
    } else if (Survey) {
      if (ThisTask == NTask-1) {
        if(ixneigh == nx) {
          ixneigh = nx-1;
          ix = nx-2;
        } 
      } else { 
        if(ixneigh == nx+1) {
          ixneigh = nx;
          ix = nx-1;
        } 
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

    density[iz+2*(NZ/2+1)*(iy+NY*ix)]           += w*itx*ity*itz;
    density[izneigh+2*(NZ/2+1)*(iy+NY*ix)]      += w*itx*ity*idz;
    density[iz+2*(NZ/2+1)*(iyneigh+NY*ix)]      += w*itx*idy*itz;
    density[izneigh+2*(NZ/2+1)*(iyneigh+NY*ix)] += w*itx*idy*idz;

    density[iz+2*(NZ/2+1)*(iy+NY*ixneigh)]           += w*idx*ity*itz;
    density[izneigh+2*(NZ/2+1)*(iy+NY*ixneigh)]      += w*idx*ity*idz;
    density[iz+2*(NZ/2+1)*(iyneigh+NY*ixneigh)]      += w*idx*idy*itz;
    density[izneigh+2*(NZ/2+1)*(iyneigh+NY*ixneigh)] += w*idx*idy*idz;

  // Triangular-Shaped Cloud interpolation
  } else if (InterpOrder == 3) {

    if (Periodic) if (x < dx) x += XMAX-XMIN;
    if (ThisTask == NTask-1) {
      if ((x < xmin+dx) || (x >= xmax+dx)) return 0.0;
    } else { 
      if ((x < xmin+dx) || (x > xmax+dx)) return 0.0;
    }

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
      }
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
    } else if (Survey) {
      if (ThisTask == NTask-1) {
        if(ixneighhi == nx) {
          ixneighhi = nx-1;
          ix = nx-2;
          ixneighlow = nx-3;
        } 
      } else { 
        if(ixneighhi == nx+1) {
          ixneighhi = nx;
          ix = nx-1;
          ixneighlow = nx-2;
        } 
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

    density[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)] += w*isx*isy*isz;
    density[iz+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)]         += w*isx*isy*itz;
    density[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ixneighlow)]  += w*isx*isy*idz;
    density[izneighlow+2*(NZ/2+1)*(iy+NY*ixneighlow)] += w*isx*ity*isz;
    density[iz+2*(NZ/2+1)*(iy+NY*ixneighlow)]         += w*isx*ity*itz;
    density[izneighhi+2*(NZ/2+1)*(iy+NY*ixneighlow)]  += w*isx*ity*idz;
    density[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)] += w*isx*idy*isz;
    density[iz+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)]         += w*isx*idy*itz;
    density[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ixneighlow)]  += w*isx*idy*idz;

    density[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ix)] += w*itx*isy*isz;
    density[iz+2*(NZ/2+1)*(iyneighlow+NY*ix)]         += w*itx*isy*itz;
    density[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ix)]  += w*itx*isy*idz;
    density[izneighlow+2*(NZ/2+1)*(iy+NY*ix)] += w*itx*ity*isz;
    density[iz+2*(NZ/2+1)*(iy+NY*ix)]         += w*itx*ity*itz;
    density[izneighhi+2*(NZ/2+1)*(iy+NY*ix)]  += w*itx*ity*idz;
    density[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ix)] += w*itx*idy*isz;
    density[iz+2*(NZ/2+1)*(iyneighhi+NY*ix)]         += w*itx*idy*itz;
    density[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ix)]  += w*itx*idy*idz;

    density[izneighlow+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)] += w*idx*isy*isz;
    density[iz+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]         += w*idx*isy*itz;
    density[izneighhi+2*(NZ/2+1)*(iyneighlow+NY*ixneighhi)]  += w*idx*isy*idz;
    density[izneighlow+2*(NZ/2+1)*(iy+NY*ixneighhi)] += w*idx*ity*isz;
    density[iz+2*(NZ/2+1)*(iy+NY*ixneighhi)]         += w*idx*ity*itz;
    density[izneighhi+2*(NZ/2+1)*(iy+NY*ixneighhi)]  += w*idx*ity*idz;
    density[izneighlow+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)] += w*idx*idy*isz;
    density[iz+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]         += w*idx*idy*itz;
    density[izneighhi+2*(NZ/2+1)*(iyneighhi+NY*ixneighhi)]  += w*idx*idy*idz;
  }
  return w;
}