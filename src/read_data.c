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

// Read in an ASCII file containing the data. 
// Every processor reads in the file but only stores the relevant parts.
// =====================================================================
double read_data_serial_ascii(void) {

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
  // We skip over any lines (starting with `#' as we assume they are headers)
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
    if (DoInterlacing) add_to_grid(tx+dx/2.0, ty+dy/2.0, tz+dz/2.0, w, XMIN_LOCAL, XMAX_LOCAL, Local_nxtra, ddg_2);

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
      ierr = MPI_Sendrecv(&(ddg_2[last_slice]),InterpOrder*alloc_slice,MPI_DOUBLE,RightTask,0,
                          &(temp_ddg[0]),InterpOrder*alloc_slice,MPI_DOUBLE,LeftTask,0,MPI_COMM_WORLD,&status);
      for (int i=0;i<InterpOrder*alloc_slice;i++) ddg_2[i] += temp_ddg[i];
    }
    free(temp_ddg);
  }

  return NGRID;

}

// Read in multiple ASCII files containing the data. 
// Every processor reads in a unique file and the processors sort out where
// it needs to go to be allocated to a grid before transferring approriately.
// ==========================================================================
double read_data_parallel_ascii(void) {

	if (ThisTask == 0) printf("\nERROR: read_data_parallel_ascii not implemented yet.\n");
    FatalError((char *)"read_param.c", 378);

    return 0.0;

}


// Add a particle to the grid using the given order of interpolation. 
// Written in such a way we only have to create and copy buffer slices of the grid on the right hand edge.
// =======================================================================================================
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
    if(ixneigh >= nx) ixneigh = nx;
    if(iyneigh >= NY) {
      iyneigh = 0;
      if (iy >= NY) iy = NY-1;
    }
    if(izneigh >= NZ) {
      izneigh = 0;
      if (iz >= NZ) iz = NZ-1;
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

    if (x < dx) x += XMAX-XMIN;
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