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

/* ======================================================================*/
/* This file contains all the prototypes for functions used in the code. */
/* ======================================================================*/

// main.c
void create_grids(void);
void destroy_grids(void);
void compute_survey_power(void);
void compute_periodic_power(void);
void output_power(double shot, double norm);
void FatalError(char * filename, int linenum);
void assign_survey_power(int multipole, int ii, int jj, int kk);

// read_data.c
void compute_nbar(int parallel, unsigned long long NDATA, unsigned long long NRAND);
void compute_fkp(unsigned long long NOBJ, struct survey_data * inputdata);
double f(double z, void *p);
double comoving_distance(double red);
double read_periodic_serial_ascii(void);
double read_survey_serial_ascii(char *inputfile, struct survey_data * inputdata);
double assign_survey_data(unsigned long long NOBJ, struct survey_data * inputdata, double prefactor);
double add_to_grid(double x, double y, double z, double w, double xmin, double xmax, int nx, double * density);

// read_params.c
void check_inputs(void);
void read_parameterfile(char * fname);
