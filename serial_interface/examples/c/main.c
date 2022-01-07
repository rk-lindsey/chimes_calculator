/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Nir Goldman (2020) 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "wrapper-C.h"

int main (int argc, char **argv) 
{
	int small = 0;	// 0 = false, 1 = true
	
  if ((argc != 4)&&(argc !=3))
  {
    printf("To run: ./test.x <parameter file> <xyz config. file>\n");
	printf("or\n");
	printf("./test.x <parameter file> <xyz config. file> <allow_replicates(0/1)>\n");
    printf("Exiting code.\n");
    exit(0);
  }
  if(argc==4)
	  small = atoi(argv[3]);

  const double GPa = 6.9479; // convert kcal/mol.A^3 to GPa
  FILE *fconf;
  int natom, i, j, k, l;
  char junk;
  double lx, ly, lz;
  double ca[3], cb[3], cc[3];
  fconf = fopen (argv[2],"r");
  fscanf(fconf,"%d\n",&natom);
  fscanf(fconf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&ca[0],&ca[1],&ca[2],&cb[0],&cb[1],&cb[2],&cc[0],&cc[1],&cc[2]);
  double vol = lx*ly*lz;
  double *stress = (double *) calloc(9,sizeof(double)); 
  double energy = 0.0;
  double xij, yij, zij, rij, dr[3];
  double xc[natom], yc[natom], zc[natom];
  double fx[natom], fy[natom], fz[natom];
  char atom_types[natom][50], *atom[natom];
  char *atype2b[2], *atype3b[3], *atype4b[4];
   
  for (i = 0; i < natom; i++) 
  {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }
  for (i = 0; i < natom; i++) 
  {
    fscanf(fconf,"%s %lf %lf %lf\n",atom_types[i],&xc[i],&yc[i],&zc[i]);
    atom[i] = atom_types[i];
  }
  fclose(fconf);
  set_chimes(small);
  
  printf("Read args:\n");
  for (i=1; i<argc; i++)
	  printf("%i %s\n",i, argv[i]);

  int rank = 0;
  init_chimes(argv[1],&rank);
  
  calculate_chimes(natom, xc, yc, zc, atom, ca, cb, cc, &energy, fx, fy, fz, stress);
  for (i = 0; i < 9; i++) {
    stress[i] *= GPa;
  }
  
  printf("\n%s\n", "Success!");
  printf("%s\t%f\n","Energy (kcal/mol):", energy); 
  printf("%s\n",   "Stress tensors (GPa)");
  printf("%s %f\n","\ts_xx: ",stress[0]);
  printf("%s %f\n","\ts_yy: ",stress[4]);
  printf("%s %f\n","\ts_zz: ",stress[8]);
  printf("%s %f\n","\ts_xy: ",stress[1]);
  printf("%s %f\n","\ts_xz: ",stress[2]);
  printf("%s %f\n","\ts_yz: ",stress[5]);
  printf("%s\n",   "Forces (kcal/mol/A)");
  for (i = 0; i <natom; i++)
  	printf("\t%f\t%f\t%f\n", fx[i], fy[i], fz[i]);
  printf("\n");
  
  #if DEBUG==1
  
  FILE *fout;
  fout = fopen("debug.dat","w");
  fprintf(fout,"%0.6f\n", energy);
  fprintf(fout,"%0.6f\n", stress[0]);
  fprintf(fout,"%0.6f\n", stress[4]);
  fprintf(fout,"%0.6f\n", stress[8]);
  fprintf(fout,"%0.6f\n", stress[1]);
  fprintf(fout,"%0.6f\n", stress[2]);
  fprintf(fout,"%0.6f\n", stress[5]);
  
  for (i = 0; i <natom; i++)
  	fprintf(fout,"%0.6f\n%0.6e\n%0.6e\n", fx[i], fy[i],fz[i]);
  fclose(fout);
  
  
  #endif

}
