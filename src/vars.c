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

/* =======================================================================================*/
/* This file contains the initialiastion of all the external variables in the header file.*/
/* =======================================================================================*/

#include "vars.h"

// MPI variables
int ierr;             // The return value for mpi routines
int NTask;            // The total number of tasks
int ThisTask;         // The rank of each task
int LeftTask;         // The first neighbouring task on the left containing particles
int RightTask;        // The first neighbouring task on the right containing particles
MPI_Status status;    // The MPI error status

// Global variables for the grids
int * Slab_to_task;       // The task to which each slice is assigned
int * Local_nx_table;     // The number of slices on each of the tasks
double * ddg;             // The density grids
double * ddg_2;           // The density grids
ptrdiff_t Local_nx;       // The number of slices on the task
ptrdiff_t last_slice;     // The last slice of the density/force grids (maybe equal to alloc_local)
ptrdiff_t Total_size;     // The total byte-size of the grids on each processor
ptrdiff_t Local_nxtra;    // The number of slices on the task accounting for buffer required for interlacing/interpolation
ptrdiff_t alloc_slice;    // The byte-size of a slice of the density/force grids
ptrdiff_t Local_x_start;  // The global start of the slices on the task
fftw_plan plan, plan_2;   // The plans for the in-place FFT of the density grid

// Parameters for input and output files
char FileBase[500];      // The base input filename
char InputDir[500];      // The input directory
char OutputDir[500];     // The output directory
int NumFilesInParallel;  // The maximum number of files to be read-in in parallel

// Parameters for grid size, particle assignment and line-of-sight read from params file
int LOS;               // The direction of the line-of-sight
int DoInterlacing;     // Whether or not to interlace the grids to reduce aliasing.
int InterpOrder;       // The interpolation order for assigning data to grids.
unsigned long long NX; // Number of grid cells in x direction
unsigned long long NY; // Number of grid cells in y direction
unsigned long long NZ; // Number of grid cells in z direction
double dx, dy, dz;     // The grid spacing
double XMIN, XMAX;	   // Physical extent of the grid in the x direction
double YMIN, YMAX;     // Physical extent of the grid in the y direction
double ZMIN, ZMAX;     // Physical extent of the grid in the z direction

// Parameters for output binning read from params file
int NK;          // Number of output k-bins
int NMU;         // Number of output mu-bins (if output_2D) is specified in params file
int Output2D;   // Whether or not to output in separate mu bins, or just the multipoles
int OutputLog;  // Whether or not to output in log k-bins
double Mink;     // minimum k to be binned
double Maxk;     // maximum k to be binned
double binwidth; // The width of the k-bins

// Storage for the actual power spectra
int * Nmodes, * Nmodes_2D;           // The number of modes in each bin
double * Pk0, * Pk2, * Pk4, * Pk_2D; // The power spectra; multipoles and 2D