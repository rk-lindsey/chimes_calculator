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
  if (argc < 3) {
    printf("To run: ./test_cpp_and_c_wrapper <parameter file> <xyz config. file>\n");
    printf("Exiting code.\n");
    exit(0);
  }
  // stress tensor calculation in chimesFF computes r.F, but leaves out volume term
  // include 1/V at the end of the calculation
  const double GPa = 6.9479; // convert kcal/mol.A^3 to GPa
  FILE *fconf;
  int natom, i, j, k, l;
  char junk;
  double lx, ly, lz, l_dummy;
  fconf = fopen (argv[2],"r");
  fscanf(fconf,"%d\n",&natom);
  fscanf(fconf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&lx,&l_dummy,&l_dummy,&l_dummy,&ly,&l_dummy,&l_dummy,&l_dummy,&lz);
  double vol = lx*ly*lz;
  double *stress = (double *) calloc(9,sizeof(double)); 
  double epot = 0.0;
  double xij, yij, zij, rij, dr[3];
  double xc[natom], yc[natom], zc[natom];
  double f2b[2][3];
  double f3b[3][3];
  double f4b[4][3];
  double ftot[natom][3];
  char atom_type[natom][50], *atom;
  char *atype2b[2], *atype3b[3], *atype4b[4]; 
  // initiliaze f2b arrays
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 3; j++) {
      f2b[i][j] = 0.0;
    }
  }
  for (i = 0; i < natom; i++) {
    for (j = 0; j < 3; j++) {
      ftot[i][j] = 0.0;
    }
  }
  for (i = 0; i < natom; i++) { // note: all atoms are C
    fscanf(fconf,"%s %lf %lf %lf\n",atom_type[i],&xc[i],&yc[i],&zc[i]);
  }
  printf("%f, %f, %f\n",lx,ly,lz);
  
  fclose(fconf);
  set_chimes();
  init_chimes(0);
  chimes_read_params(argv[1]);
  int order2b, order3b, order4b;
  order2b = get_chimes_2b_order();
  order3b = get_chimes_3b_order();
  order4b = get_chimes_4b_order();
// compute 2B interaction 
  double rcut_2b;
  rcut_2b = get_chimes_max_2b_cutoff();
  for (i = 0; i < natom-1; i++) {
    for (j = i+1; j < natom; j++) {
     // compute relative coordinates and apply minimum image PBC
      xij = (xc[i] - xc[j]);
      xij -= lx*round(xij/lx);
      dr[0] = xij;
      yij = (yc[i] - yc[j]);
      yij -= ly*round(yij/ly);
      dr[1] = yij;
      zij = (zc[i] - zc[j]);
      zij -= lz*round(zij/lz);
      dr[2] = zij;
      rij = sqrt(xij*xij + yij*yij + zij*zij);
      f2b[0][0] = ftot[i][0];
      f2b[0][1] = ftot[i][1];
      f2b[0][2] = ftot[i][2];
      f2b[1][0] = ftot[j][0];
      f2b[1][1] = ftot[j][1];
      f2b[1][2] = ftot[j][2];
      atype2b[0] = atom_type[i];
      atype2b[1] = atom_type[j];
      if (rij < rcut_2b) {
        // f2b, stress tensor, epot are all cumulative
        chimes_compute_2b_props(rij, dr, atype2b, f2b, stress, &epot);
      }
      ftot[i][0] = f2b[0][0];
      ftot[i][1] = f2b[0][1];
      ftot[i][2] = f2b[0][2];
      ftot[j][0] = f2b[1][0];
      ftot[j][1] = f2b[1][1];
      ftot[j][2] = f2b[1][2];
    }
  }
  printf("2B RESULTS....\n");
  printf("epot = %lf\n",epot);
// compute 3B interaction 
  double rcut_3b;
  rcut_3b = get_chimes_max_3b_cutoff();
  if (order3b > 0) 
  {
  // initiliaze f3b array
    for (i = 0; i < 3; i++) 
    {
      for (j = 0; j < 3; j++) 
      {
        f3b[i][j] = 0.0;
      }
    }
    double dist_3b[3][3];
    double dr_3b[3];
    double rik,rjk;
    double xik,yik,zik;
    double xjk,yjk,zjk;
    for (i = 0; i < natom-2; i++) {
      for (j = i+1; j < natom-1; j++) {
        for (k = j+1; k < natom; k++) {
       // compute relative coordinates and apply minimum image PBC
       // order in chimesFF is ij, ik, jk
       //ij pairs 
          xij = (xc[i] - xc[j]);
          xij -= lx*round(xij/lx);
          yij = (yc[i] - yc[j]);
          yij -= ly*round(yij/ly);
          zij = (zc[i] - zc[j]);
          zij -= lz*round(zij/lz);
          dist_3b[0][0] = xij;
          dist_3b[0][1] = yij;
          dist_3b[0][2] = zij;
          dr_3b[0] = sqrt(xij*xij + yij*yij + zij*zij);
       //ik pairs
          xik = (xc[i] - xc[k]);
          xik -= lx*round(xik/lx);
          yik = (yc[i] - yc[k]);
          yik -= ly*round(yik/ly);
          zik = (zc[i] - zc[k]);
          zik -= lz*round(zik/lz);
          dist_3b[1][0] = xik;
          dist_3b[1][1] = yik;
          dist_3b[1][2] = zik;
          dr_3b[1] = sqrt(xik*xik + yik*yik + zik*zik);
       //jk pairs 
          xjk = (xc[j] - xc[k]);
          xjk -= lx*round(xjk/lx);
          yjk = (yc[j] - yc[k]);
          yjk -= ly*round(yjk/ly);
          zjk = (zc[j] - zc[k]);
          zjk -= lz*round(zjk/lz);
          dist_3b[2][0] = xjk;
          dist_3b[2][1] = yjk;
          dist_3b[2][2] = zjk;
          dr_3b[2] = sqrt(xjk*xjk + yjk*yjk + zjk*zjk);
          f3b[0][0] = ftot[i][0];
          f3b[0][1] = ftot[i][1];
          f3b[0][2] = ftot[i][2];
          f3b[1][0] = ftot[j][0];
          f3b[1][1] = ftot[j][1];
          f3b[1][2] = ftot[j][2];
          f3b[2][0] = ftot[k][0];
          f3b[2][1] = ftot[k][1];
          f3b[2][2] = ftot[k][2];
          // array of atom types for this triplet
          atype3b[0] = atom_type[i];
          atype3b[1] = atom_type[j];
          atype3b[2] = atom_type[k];
          if (dr_3b[0] < rcut_3b) {
            if (dr_3b[1] < rcut_3b) {
              if (dr_3b[2] < rcut_3b) {
                chimes_compute_3b_props(dr_3b, dist_3b, atype3b, f3b, stress, &epot);
              }
            }
          }
          ftot[i][0] = f3b[0][0];
          ftot[i][1] = f3b[0][1];
          ftot[i][2] = f3b[0][2];
          ftot[j][0] = f3b[1][0];
          ftot[j][1] = f3b[1][1];
          ftot[j][2] = f3b[1][2];
          ftot[k][0] = f3b[2][0];
          ftot[k][1] = f3b[2][1];
          ftot[k][2] = f3b[2][2];
        }
      }
    }
    printf("+3B RESULTS....\n");
    printf("epot = %lf\n",epot);
  }
// compute 4B interaction 
  double rcut_4b;
  rcut_4b = get_chimes_max_4b_cutoff();
  if (order4b > 0) {
  // initiliaze f4b array
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 3; j++) {
        f4b[i][j] = 0.0;
      }
    }
    double dist_4b[6][3];
    double dr_4b[6];
    double xik,yik,zik;
    double xil,yil,zil;
    double xjk,yjk,zjk;
    double xjl,yjl,zjl;
    double xkl,ykl,zkl;
    for (i = 0; i < natom-3; i++) {
      for (j = i+1; j < natom-2; j++) {
        for (k = j+1; k < natom-1; k++) {
          for (l = k+1; l < natom; l++) {
         // compute relative coordinates and apply minimum image PBC
         // order in chimesFF is: ij, ik, il, jk, jl, kl
         //ij pairs 
            xij = (xc[i] - xc[j]);
            xij -= lx*round(xij/lx);
            yij = (yc[i] - yc[j]);
            yij -= ly*round(yij/ly);
            zij = (zc[i] - zc[j]);
            zij -= lz*round(zij/lz);
            dist_4b[0][0] = xij;
            dist_4b[0][1] = yij;
            dist_4b[0][2] = zij;
            dr_4b[0] = sqrt(xij*xij + yij*yij + zij*zij);
         //ik pairs 
            xik = (xc[i] - xc[k]);
            xik -= lx*round(xik/lx);
            yik = (yc[i] - yc[k]);
            yik -= ly*round(yik/ly);
            zik = (zc[i] - zc[k]);
            zik -= lz*round(zik/lz);
            dist_4b[1][0] = xik;
            dist_4b[1][1] = yik;
            dist_4b[1][2] = zik;
            dr_4b[1] = sqrt(xik*xik + yik*yik + zik*zik);
         //il pairs 
            xil = (xc[i] - xc[l]);
            xil -= lx*round(xil/lx);
            yil = (yc[i] - yc[l]);
            yil -= ly*round(yil/ly);
            zil = (zc[i] - zc[l]);
            zil -= lz*round(zil/lz);
            dist_4b[2][0] = xil;
            dist_4b[2][1] = yil;
            dist_4b[2][2] = zil;
            dr_4b[2] = sqrt(xil*xil + yil*yil + zil*zil);
         //jk pairs 
            xjk = (xc[j] - xc[k]);
            xjk -= lx*round(xjk/lx);
            yjk = (yc[j] - yc[k]);
            yjk -= ly*round(yjk/ly);
            zjk = (zc[j] - zc[k]);
            zjk -= lz*round(zjk/lz);
            dist_4b[3][0] = xjk;
            dist_4b[3][1] = yjk;
            dist_4b[3][2] = zjk;
            dr_4b[3] = sqrt(xjk*xjk + yjk*yjk + zjk*zjk);
         //jl pairs 
            xjl = (xc[j] - xc[l]);
            xjl -= lx*round(xjl/lx);
            yjl = (yc[j] - yc[l]);
            yjl -= ly*round(yjl/ly);
            zjl = (zc[j] - zc[l]);
            zjl -= lz*round(zjl/lz);
            dist_4b[4][0] = xjl;
            dist_4b[4][1] = yjl;
            dist_4b[4][2] = zjl;
            dr_4b[4] = sqrt(xjl*xjl + yjl*yjl + zjl*zjl);
         //kl pairs 
            xkl = (xc[k] - xc[l]);
            xkl -= lx*round(xkl/lx);
            ykl = (yc[k] - yc[l]);
            ykl -= ly*round(ykl/ly);
            zkl = (zc[k] - zc[l]);
            zkl -= lz*round(zkl/lz);
            dist_4b[5][0] = xkl;
            dist_4b[5][1] = ykl;
            dist_4b[5][2] = zkl;
            dr_4b[5] = sqrt(xkl*xkl + ykl*ykl + zkl*zkl);
          // set up f4b array
            f4b[0][0] = ftot[i][0];
            f4b[0][1] = ftot[i][1];
            f4b[0][2] = ftot[i][2];
            f4b[1][0] = ftot[j][0];
            f4b[1][1] = ftot[j][1];
            f4b[1][2] = ftot[j][2];
            f4b[2][0] = ftot[k][0];
            f4b[2][1] = ftot[k][1];
            f4b[2][2] = ftot[k][2];
            f4b[3][0] = ftot[l][0];
            f4b[3][1] = ftot[l][1];
            f4b[3][2] = ftot[l][2];
            atype4b[0] = atom_type[i];
            atype4b[1] = atom_type[j];
            atype4b[2] = atom_type[k];
            atype4b[3] = atom_type[l];
          // perform chimes 4B calc.
            if (dr_4b[0] < rcut_4b) {
              if (dr_4b[1] < rcut_4b) {
                if (dr_4b[2] < rcut_4b) {
                  if (dr_4b[3] < rcut_4b) {
                    if (dr_4b[4] < rcut_4b) {
                      if (dr_4b[5] < rcut_4b) {
                        chimes_compute_4b_props(dr_4b, dist_4b, atype4b, f4b, stress, &epot);
                      }
                    }
                  }
                }
              }
            }
            ftot[i][0] = f4b[0][0];
            ftot[i][1] = f4b[0][1];
            ftot[i][2] = f4b[0][2];
            ftot[j][0] = f4b[1][0];
            ftot[j][1] = f4b[1][1];
            ftot[j][2] = f4b[1][2];
            ftot[k][0] = f4b[2][0];
            ftot[k][1] = f4b[2][1];
            ftot[k][2] = f4b[2][2];
            ftot[l][0] = f4b[3][0];
            ftot[l][1] = f4b[3][1];
            ftot[l][2] = f4b[3][2];
          }
        }
      }
    }
    printf("+4B RESULTS....\n");
    printf("epot = %lf\n",epot);
  }
  // divide stess by 1/V and convert to GPa
  for (i = 0; i < 9; i++) {
    stress[i] /= vol;
    stress[i] *= GPa;
  }
  printf("Printing final forces in output_lib.xyzf\n");
  FILE *fout;
  fout = fopen("output_lib.xyzf","w");
  fprintf(fout,"%d\n",natom);
  fprintf(fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",lx,ly,lz,stress[0],stress[1],stress[2],stress[3],stress[4],stress[5],stress[6],stress[7],stress[8],epot);
  for (i = 0; i < natom; i++) {
    fprintf(fout,"%s %lf %lf %lf %lf %lf %lf\n",atom_type[i],xc[i],yc[i],zc[i],ftot[i][0],ftot[i][1],ftot[i][2]);
  }
}
