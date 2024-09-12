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

/* =========================================================================*/
/* This file contains routines to read in and check the run parameters file.*/
/* =========================================================================*/

#include "vars.h"
#include "proto.h"

// Read in the list of input parameters
// ====================================
void read_parameterfile(char * fname) {

#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  char buf[500],buf1[500],buf2[500],buf3[500];
  int i,j,nt;
  int id[MAXTAGS];
  int errorFlag = 0;

  // read parameter file on all processes for simplicity

  nt = 0;

  //strcpy(tag[nt], "InputDir");
  //addr[nt] = InputDir;
  //id[nt++] = STRING;

  //strcpy(tag[nt], "OutputDir");
  //addr[nt] = OutputDir;
  //id[nt++] = STRING;

  //strcpy(tag[nt], "FileBase");
  //addr[nt] = FileBase;
  //id[nt++] = STRING;

  //strcpy(tag[nt], "RandFileBase");
  //addr[nt] = RandFileBase;
  //id[nt++] = STRING;

  strcpy(tag[nt], "NumFilesInParallel");
  addr[nt] = &NumFilesInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "NX");
  addr[nt] = &NX;
  id[nt++] = INT;

  strcpy(tag[nt], "NY");
  addr[nt] = &NY;
  id[nt++] = INT;

  strcpy(tag[nt], "NZ");
  addr[nt] = &NZ;
  id[nt++] = INT;

  strcpy(tag[nt], "XMIN");
  addr[nt] = &XMIN;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "XMAX");
  addr[nt] = &XMAX;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "YMIN");
  addr[nt] = &YMIN;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "YMAX");
  addr[nt] = &YMAX;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ZMIN");
  addr[nt] = &ZMIN;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ZMAX");
  addr[nt] = &ZMAX;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "DoInterlacing");
  addr[nt] = &DoInterlacing;
  id[nt++] = INT;

  strcpy(tag[nt], "InterpOrder");
  addr[nt] = &InterpOrder;
  id[nt++] = INT;

  strcpy(tag[nt], "NK");
  addr[nt] = &NK;
  id[nt++] = INT;

  strcpy(tag[nt], "NMU");
  addr[nt] = &NMU;
  id[nt++] = INT;

  strcpy(tag[nt], "Output2D");
  addr[nt] = &Output2D;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputLog");
  addr[nt] = &OutputLog;
  id[nt++] = INT;

  strcpy(tag[nt], "Mink");
  addr[nt] = &Mink;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Maxk");
  addr[nt] = &Maxk;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omega_m");
  addr[nt] = &Omega_m;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Redshift");
  addr[nt] = &Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Periodic");
  addr[nt] = &Periodic;
  id[nt++] = INT;

  strcpy(tag[nt], "LOS");
  addr[nt] = &LOS;
  id[nt++] = INT;

  strcpy(tag[nt], "Survey");
  addr[nt] = &Survey;
  id[nt++] = INT;

  strcpy(tag[nt], "Momentum");
  addr[nt] = &Momentum;
  id[nt++] = INT;

  strcpy(tag[nt], "X_Origin");
  addr[nt] = &X_Origin;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Y_Origin");
  addr[nt] = &Y_Origin;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Z_Origin");
  addr[nt] = &Z_Origin;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "REDMININ");
  addr[nt] = &REDMININ;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "REDMAXIN");
  addr[nt] = &REDMAXIN;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "NOBJ_Max");
  addr[nt] = &NOBJ_Max;
  id[nt++] = INT;

  strcpy(tag[nt], "Coord_Type");
  addr[nt] = &Coord_Type;
  id[nt++] = INT;

  strcpy(tag[nt], "x_Column");
  addr[nt] = &x_Column;
  id[nt++] = INT;

  strcpy(tag[nt], "y_Column");
  addr[nt] = &y_Column;
  id[nt++] = INT;

  strcpy(tag[nt], "z_Column");
  addr[nt] = &z_Column;
  id[nt++] = INT;

  strcpy(tag[nt], "NBAR_Column");
  addr[nt] = &NBAR_Column;
  id[nt++] = INT;

  strcpy(tag[nt], "SkyArea");
  addr[nt] = &SkyArea;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "FKP_Column");
  addr[nt] = &FKP_Column;
  id[nt++] = INT;

  strcpy(tag[nt], "FKP_Pk");
  addr[nt] = &FKP_Pk;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "FKP_Pk_mom");
  addr[nt] = &FKP_Pk_mom;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "FKP_Pk_cross");
  addr[nt] = &FKP_Pk_cross;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Odd_Multipoles");
  addr[nt] = &Odd_Multipoles;
  id[nt++] = INT;

  if((fd = fopen(fname, "r"))) {
    fflush(stdout);
    while(!feof(fd)) {
      buf[0] = 0;
      fgets(buf, 500, fd);

      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
      if(buf1[0] == '%') continue;

      for(i = 0, j = -1; i < nt; i++) {
        if(strcmp(buf1, tag[i]) == 0)  {
          j = i;
	        tag[i][0] = 0;
	        break;
	      }
      }
      
      if(j >= 0) {
  	    switch (id[j]) {
  	      case FLOAT:
  	        *((double *) addr[j]) = atof(buf2);
  	        break;
  	      case STRING:
  	        strcpy((char *)addr[j], buf2);
  	        break;
  	      case INT:
  	        *((int *) addr[j]) = atoi(buf2);
  	        break;
  	    }
      } else {
        if(ThisTask == 0) fprintf(stdout,"\nERROR: In file %s:  Tag '%s' not allowed or multiple defined.\n\n",fname, buf1);
	      errorFlag = 1;
      }
    }
    fclose(fd);
    for(i = 0; i < nt; i++) {
      if(*tag[i]) {
        if(ThisTask == 0) fprintf(stdout, "\nERROR: I miss a value for tag '%s' in parameter file '%s'.\n\n", tag[i], fname);
        errorFlag = 1;
      }
    }
  } else {
    if(ThisTask == 0) fprintf(stdout,"\nERROR: Parameter file '%s' not found.\n\n",fname);
    FatalError((char *)"read_param.c", 330);
  }

  if(errorFlag) {
    fflush(stdout);
    MPI_Finalize();
    exit(1);
  }

  check_inputs();

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS

  return;

}

// Check the run parameters to ensure suitable values are used.
// ============================================================
void check_inputs(void) {

  if (NumFilesInParallel < 1) {
    if (ThisTask == 0) {
      printf("\nERROR: `NumFilesInParallel' is less than 1 so no processors will be writing out the data.\n");
      printf("        Please set NumFileInParallel > 0'.\n\n");
    }
    FatalError((char *)"read_param.c", 297);
  }

  if (NTask < NumFilesInParallel) {
    if (ThisTask == 0) {
      printf("\nWARNING: Number of processors smaller than `NumFilesInParallel'.\n");
      printf("         Setting NumFilesInParallel = Number of processors'.\n\n");
    }
  }

  if (Periodic) {
    if ((LOS != 0) && (LOS != 1) && (LOS != 2) && (LOS != 3)) {
      if (ThisTask == 0) {
        printf("\nERROR: `LOS' is neither 0 (real-space), 1, 2 nor 3.\n");
        printf("        Please choose a suitable value for `LOS'.\n\n");
      }
      FatalError((char *)"read_param.c", 313);
    }
  }

  if ((InterpOrder != 1) && (InterpOrder != 2) && (InterpOrder != 3)) {
    if (ThisTask == 0) {
      printf("\nERROR: `InterpOrder' is neither 1, 2 nor 3.\n");
      printf("        Please choose a suitable value for `InterpOrder'.\n\n");
    }
    FatalError((char *)"read_param.c", 329);
  }

  if ((DoInterlacing != 0) && (DoInterlacing != 1)) {
    if (ThisTask == 0) {
      printf("\nERROR: `DoInterlacing' is neither 0 nor 1.\n");
      printf("        Please choose a suitable value for `DoInterlacing'.\n\n");
    }
    FatalError((char *)"read_param.c", 337);
  }

  if ((Output2D != 0) && (Output2D != 1)) {
    if (ThisTask == 0) {
      printf("\nERROR: `Output2D' is neither 0 nor 1.\n");
      printf("        Please choose a suitable value for `Output2D'.\n\n");
    }
    FatalError((char *)"read_param.c", 345);
  }

  if (Output2D) {
    if (NMU <= 0) {
      if (ThisTask == 0) {
        printf("\nERROR: `NMU' is less than 1.  if Output2D is 1, we must choose some mu-bins to output.\n");
        printf("        Please set 'NMU' > 0.\n\n");
      }
      FatalError((char *)"read_param.c", 354);
    }
  }

  if ((OutputLog != 0) && (OutputLog != 1)) {
    if (ThisTask == 0) {
      printf("\nERROR: `OutputLog' is neither 0 nor 1.\n");
      printf("        Please choose a suitable value for `OutputLog'.\n\n");
    }
    FatalError((char *)"read_param.c", 363);
  }

  if (NX <= 0) {
    if (ThisTask == 0) {
      printf("\nERROR: `NX' is less than 1. We must generate a mesh.\n");
      printf("        Please set 'NX' > 0.\n\n");
    }
    FatalError((char *)"read_param.c", 371);
  }

  if (NY <= 0) {
    if (ThisTask == 0) {
      printf("\nERROR: `NY' is less than 1. We must generate a mesh.\n");
      printf("        Please set 'NY' > 0.\n\n");
    }
    FatalError((char *)"read_param.c", 379);
  }

  if (NZ <= 0) {
    if (ThisTask == 0) {
      printf("\nERROR: `NZ' is less than 1. We must generate a mesh.\n");
      printf("        Please set 'NZ' > 0.\n\n");
    }
    FatalError((char *)"read_param.c", 387);
  }

  if (NK <= 0) {
    if (ThisTask == 0) {
      printf("\nERROR: `NK' is less than 1. We must choose some k-bins to output.\n");
      printf("        Please set 'NK' > 0.\n\n");
    }
    FatalError((char *)"read_param.c", 395);
  }

  if (XMAX < XMIN) {
    if (ThisTask == 0) {
      printf("\nERROR: `XMAX' is less than `XMIN'.\n");
      printf("        Please set `XMAX' > `XMIN'.\n\n");
    }
    FatalError((char *)"read_param.c", 403);
  }

  if (YMAX < YMIN) {
    if (ThisTask == 0) {
      printf("\nERROR: `XMAX' is less than `XMIN'.\n");
      printf("        Please set `XMAX' > `XMIN'.\n\n");
    }
    FatalError((char *)"read_param.c", 411);
  }

  if (ZMAX < ZMIN) {
    if (ThisTask == 0) {
      printf("\nERROR: `XMAX' is less than `XMIN'.\n");
      printf("        Please set `XMAX' > `XMIN'.\n\n");
    }
    FatalError((char *)"read_param.c", 419);
  }

  if (Maxk < Mink) {
    if (ThisTask == 0) {
      printf("\nERROR: `Maxk' is less than `Mink'.\n");
      printf("        Please set `Maxk' > `Mink'.\n\n");
    }
    FatalError((char *)"read_param.c", 427);
  }

  return;
}