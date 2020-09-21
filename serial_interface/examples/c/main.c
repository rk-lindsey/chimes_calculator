#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "wrapper-C.h"

int main (int argc, char **argv) 
{
  if (argc < 3) {
    printf("To run: ./test.x <parameter file> <xyz config. file>\n");
    printf("Exiting code.\n");
    exit(0);
  }
  const double GPa = 6.9479; // convert kcal/mol.A^3 to GPa
  FILE *fconf;
  int natom, i, j, k, l;
  int nlayer;
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
  for (i = 0; i < natom; i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }
  for (i = 0; i < natom; i++) { // note: all atoms are C
    fscanf(fconf,"%s %lf %lf %lf\n",atom_types[i],&xc[i],&yc[i],&zc[i]);
    atom[i] = atom_types[i];
  }
  fclose(fconf);
  set_chimes();
  nlayer = 4;
  init_chimes(argv[1], &nlayer);
  calculate_chimes(natom, xc, yc, zc, atom, ca, cb, cc, &energy, fx, fy, fz, stress);
  printf("ener: %lf\n",energy);
  for (i = 0; i < 9; i++) {
    stress[i] *= GPa;
  }
  printf("Printing final forces in output_lib.xyzf\n");
  FILE *fout;
  fout = fopen("output_lib_simple.xyzf","w");
  fprintf(fout,"%d\n",natom);
  fprintf(fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",ca[0],cb[1],cc[2],stress[0],stress[1],stress[2],stress[3],stress[4],stress[5],stress[6],stress[7],stress[8],energy);
  for (i = 0; i < natom; i++) {
    fprintf(fout,"%s %lf %lf %lf %lf %lf %lf\n",atom_types[i],xc[i],yc[i],zc[i],fx[i],fy[i],fz[i]);
  }
}
