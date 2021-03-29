/*
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Nir Goldman (2020)
*/

#include<vector>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<map>

using namespace std;

#include "serial_chimes_interface.h"
#include "wrapper-C.h"
static  serial_chimes_interface chimes, *chimes_ptr;

void set_chimes () {
  chimes_ptr = &chimes;
}

void init_chimes (char *param_file, int *rank) {
  chimes_ptr->init_chimesFF(param_file, *rank);
}

void calculate_chimes(int natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9])
{
  vector<double>    x_vec(natom);
  vector<double>    y_vec(natom);
  vector<double>    z_vec(natom);

  vector<vector<double> > force_vec;
  force_vec.resize(natom, vector<double>(3,0.0));


  vector<string> atom_types_vec;
  atom_types_vec.resize(natom);


  for (int i = 0; i < natom; i++) {
    x_vec[i] = xc[i];
    y_vec[i] = yc[i];
    z_vec[i] = zc[i];
    force_vec[i][0] = fx[i];
    force_vec[i][1] = fy[i];
    force_vec[i][2] = fz[i];
    atom_types_vec[i] = atom_types[i];
  }
  vector<double> stress_vec(9,0.0);
  for (int i = 0; i < 9; i++) {
    stress_vec[i] = stress[i];
  }
  //for (int i = 0; i < 9; i++) {
  //  stress[i] = stress_vec[i];
  //}
  vector<double>cell_a_vec(3);
  vector<double>cell_b_vec(3);
  vector<double>cell_c_vec(3);
  cell_a_vec[0] = ca[0];
  cell_a_vec[1] = ca[1];
  cell_a_vec[2] = ca[2];
  cell_b_vec[0] = cb[0];
  cell_b_vec[1] = cb[1];
  cell_b_vec[2] = cb[2];
  cell_c_vec[0] = cc[0];
  cell_c_vec[1] = cc[1];
  cell_c_vec[2] = cc[2];

  chimes_ptr->calculate(x_vec, y_vec, z_vec, cell_a_vec, cell_b_vec, cell_c_vec, atom_types_vec, *energy, force_vec, stress_vec);
  for (int i = 0; i < natom; i++) {
    fx[i] = force_vec[i][0];
    fy[i] = force_vec[i][1];
    fz[i] = force_vec[i][2];
  }
  for (int i = 0; i < 9; i++) {
    stress[i] = stress_vec[i];
  }
}
