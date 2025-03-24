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

/* =======================================================*/
/* This file contains all the global variable definitions.*/
/* =======================================================*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

// MPI and FFTW libraries
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

// GSL Libraries
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

// MPI variables
extern int ierr;             // The return value for mpi routines
extern int NTask;            // The total number of tasks
extern int ThisTask;         // The rank of each task
extern int LeftTask;         // The first neighbouring task on the left containing particles
extern int RightTask;        // The first neighbouring task on the right containing particles
extern MPI_Status status;    // The MPI error status

// Global variables for the grids
extern int * Slab_to_task;         // The task to which each slice is assigned
extern int * Local_nx_table;       // The number of slices on each of the tasks
extern double * F0, * F0_mom;                // Density grid for storing fourier transformed monopole
extern double * ddg, * ddg_2;      // The density grids for storing data
extern double * ddg_interlace;     // Extra density grids for interlacing
extern double * ddg_interlace_2;   // Extra density grids for interlacing
extern double * ddg_mom, * ddg_mom_2;      // The density grids for storing data
extern double * ddg_mom_interlace;     // Extra density grids for interlacing
extern double * ddg_mom_interlace_2;   // Extra density grids for interlacing
extern ptrdiff_t Local_nx;         // The number of slices on the task
extern ptrdiff_t last_slice;       // The last slice of the density/force grids (maybe equal to alloc_local)
extern ptrdiff_t Total_size;       // The total byte-size of the grids on each processor
extern ptrdiff_t Local_nxtra;      // The number of slices on the task accounting for buffer required for interlacing/interpolation
extern ptrdiff_t alloc_slice;      // The byte-size of a slice of the density/force grids
extern ptrdiff_t Local_x_start;    // The global start of the slices on the task
extern fftw_plan plan, plan_mom;     // The plans for the in-place FFT of the density grid
extern fftw_plan plan_interlace;   // The plans for the in-place FFT of the density grid
extern fftw_plan plan_mom_interlace; // The plans for the in-place FFT of the density grid

// Parameters for input and output files
extern char FileBase[500];      // The base input filename
extern char InputDir[500];      // The input directory
extern char OutputDir[500];     // The output directory
extern char RandFileBase[500];  // The random base input filename
extern int NumFilesInParallel;  // The maximum number of files to be read-in in parallel

// Parameters for grid size, particle assignment and line-of-sight read from params file
extern int DoInterlacing;     // Whether or not to interlace the grids to reduce aliasing.
extern int InterpOrder;       // The interpolation order for assigning data to grids.
extern unsigned long long NX; // Number of grid cells in x direction
extern unsigned long long NY; // Number of grid cells in y direction
extern unsigned long long NZ; // Number of grid cells in z direction
extern double dx, dy, dz;     // The grid spacing
extern double XMIN, XMAX;	  // Physical extent of the grid in the x direction
extern double YMIN, YMAX;     // Physical extent of the grid in the y direction
extern double ZMIN, ZMAX;     // Physical extent of the grid in the z direction

// Parameters for output binning read from params file
extern int NK;          // Number of output k-bins
extern int NMU;         // Number of output mu-bins (if output_2D) is specified in params file
extern int Output2D;    // Whether or not to output in separate mu bins, or just the multipoles
extern int OutputLog;   // Whether or not to output in log k-bins
extern double Mink;     // Minimum k to be binned
extern double Maxk;     // Maximum k to be binned
extern double binwidth; // The width of the k-bins

// Storage for the actual power spectra
extern int * Nmodes, * Nmodes_2D;                        // The number of modes in each bin
extern double * Pk0, *Pk1, * Pk2, * Pk3, * Pk4, * Pk_2D; // The power spectra; multipoles and 2D

// Cosmology and other constants
#define  LightSpeed  299792.458                  // The speed of light in km/s
extern double Omega_m;                           // The matter density at the present day
extern double Redshift;                          // The redshift (for periodic boxes)
extern gsl_spline * red_spline, * dist_spline;   // Splines for the redshift-distance relation
extern gsl_interp_accel * red_acc, * dist_acc;   // Splines for the redshift-distance relation

// Parameters for simulation boxes
extern int LOS;      // The direction of the line-of-sight
extern int Periodic; // Tell the code this is data from a periodic simulation (this or Survey must be set to 1)

// Parameters for survey data
extern int Survey;         // Tell the code this is data from a periodic simulation (this or Periodic must be set to 1)
extern int Momentum;       // Whether or not to compute the momentum auto/cross power spectra (0 == density, 1 == momentum, otherwise == cross)
extern int NOBJ_Max;       // The maximum number of objects to read in (of any single type type, data or randoms, for allocating space)
extern int x_Column;       // The column in the input file containing the x coordinates or RA values
extern int y_Column;       // The column in the input file containing the y coordinates or Dec values
extern int z_Column;       // The column in the input file containing the z coordinates or redshift values
extern int Coord_Type;     // Are the input coordinates in cartesian coordinates, angular degrees, or angular radians? (0, 1 and othewise respectively)
extern int FKP_Column;     // The column in the input file containing the FKP weights. Set to -1 to compute these internally using FKP_Pk value
extern int NBAR_Column;    // The column in the input file containing the nbar values. Set to -1 to compute these internally using SkyArea value
extern int Odd_Multipoles; // Flag to tell the code whether or not to compute the odd-order multipoles
extern double SkyArea;    // The effective solid angle area of the survey in square degrees (if necessary)
extern double FKP_Pk;     // The value to use for computing the FKP weights (if necessary)
extern double FKP_Pk_mom;     // The value to use for computing the FKP weights (if necessary)
extern double FKP_Pk_cross;     // The value to use for computing the FKP weights (if necessary)
extern double REDMIN;     // The minimum redshift of the data (used for computing the nbar)    
extern double REDMAX;     // The maximum redshift of the data (used for computing the nbar) 
extern double REDMININ;       // The minimum redshift to include in the input
extern double REDMAXIN;       // The maximum redshift to include in the input
extern double X_Origin;   // The coodinate of the observer in the box in the x-direction
extern double Y_Origin;   // The coodinate of the observer in the box in the y-direction
extern double Z_Origin;   // The coodinate of the observer in the box in the z-direction
extern double nsq, vr_ave, vr2_ave, vr3_ave, vr4_ave; // Variables needed for momentum power spectra 

// Structures for survey data
extern struct survey_data { 
  float coord[3];             // The cartesian coordinates of the object
  float redshift;             // The redshift of the object
  float nbar;                 // The nbar at the objects location
  float weight, weight_pv;            // The weight of the object
  float pv, pverr;        // The peculiar velocity and its error
} *data, *randoms;
