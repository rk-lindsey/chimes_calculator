/* Code author: Nir Goldman (2020) */

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

#include "chimesFF.h"
#include "wrapper-C.h"
static chimesFF chimes_start, *chimes_ptr;

double get_chimes_max_2b_cutoff() {
  double rcut_2b = chimes_ptr->max_cutoff_2B();
  return rcut_2b;
}
  
double get_chimes_max_3b_cutoff() {
  double rcut_3b = chimes_ptr->max_cutoff_3B();
  return rcut_3b;
}
  
double get_chimes_max_4b_cutoff() {
  double rcut_4b = chimes_ptr->max_cutoff_4B();
  return rcut_4b;
}

void set_chimes () {
  chimes_ptr = &chimes_start;
}

int get_chimes_2b_order() {
  int chimes2b_order = chimes_ptr->poly_orders[0];
  return chimes2b_order;
}

int get_chimes_3b_order() {
  int chimes3b_order;
  if (chimes_ptr->poly_orders.size() > 1) {
    chimes3b_order = chimes_ptr->poly_orders[1];
  } else {
    chimes3b_order = 0;
  }
  return chimes3b_order;
}

int get_chimes_4b_order() {
  int chimes4b_order;
  if (chimes_ptr->poly_orders.size() > 2) {
    chimes4b_order = chimes_ptr->poly_orders[2];
  } else {
    chimes4b_order = 0;
  }
  return chimes4b_order;
}

void init_chimes (int rank) {
  chimes_ptr->init(rank);
}

void chimes_read_params(char *param_file) {
  chimes_ptr->read_parameters(param_file);
}

void chimes_compute_2b_props(double rij, double dr[3], char *atype2b[2], double force[2][3], double stress[9], double *epot) {
  // convert all doubles, etc., from C to type vector for C++
  // declare needed vectors for chimes
  vector <double> dr_vec(3);
  dr_vec[0] = dr[0];
  dr_vec[1] = dr[1];
  dr_vec[2] = dr[2];
  vector <int> type_vec(2);
  type_vec[0] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype2b[0]));
  type_vec[1] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype2b[1]));
  //type_vec[0] = chimes_ptr->atmtoidx[atype2b[0]];
  //type_vec[1] = chimes_ptr->atmtoidx[atype2b[1]];
  vector<vector<double*> >force_vec;
  force_vec.resize(2,vector<double*>(3));
  force_vec[0][0] = &force[0][0];
  force_vec[0][1] = &force[0][1];
  force_vec[0][2] = &force[0][2];
  force_vec[1][0] = &force[1][0];
  force_vec[1][1] = &force[1][1];
  force_vec[1][2] = &force[1][2];
  vector<double*> stress_vec(9);
  stress_vec[0] = &stress[0];
  stress_vec[1] = &stress[1];
  stress_vec[2] = &stress[2];
  stress_vec[3] = &stress[3];
  stress_vec[4] = &stress[4];
  stress_vec[5] = &stress[5];
  stress_vec[6] = &stress[6];
  stress_vec[7] = &stress[7];
  stress_vec[8] = &stress[8];
  chimes_ptr->compute_2B(rij, dr_vec, type_vec, force_vec, stress_vec, *epot);
  // save forces and stress tensor components 
  force[0][0] = *force_vec[0][0];
  force[0][1] = *force_vec[0][1];
  force[0][2] = *force_vec[0][2];
  force[1][0] = *force_vec[1][0];
  force[1][1] = *force_vec[1][1];
  force[1][2] = *force_vec[1][2];
  stress[0] = *stress_vec[0];
  stress[1] = *stress_vec[1];
  stress[2] = *stress_vec[2];
  stress[3] = *stress_vec[3];
  stress[4] = *stress_vec[4];
  stress[5] = *stress_vec[5];
  stress[6] = *stress_vec[6];
  stress[7] = *stress_vec[7];
  stress[8] = *stress_vec[8];
}

void chimes_compute_3b_props(double dr_3b[3], double dist_3b[3][3], char *atype3b[3], double f3b[3][3], double stress[9], double *epot) {
  // convert all doubles, etc., from C to type vector for C++
  vector <double> dr_3b_vec(3);
  dr_3b_vec[0] = dr_3b[0];
  dr_3b_vec[1] = dr_3b[1];
  dr_3b_vec[2] = dr_3b[2];
  vector< vector<double> > dist_3b_vec(3, vector<double> (3));
  dist_3b_vec[0][0] = dist_3b[0][0];
  dist_3b_vec[0][1] = dist_3b[0][1];
  dist_3b_vec[0][2] = dist_3b[0][2];
  dist_3b_vec[1][0] = dist_3b[1][0];
  dist_3b_vec[1][1] = dist_3b[1][1];
  dist_3b_vec[1][2] = dist_3b[1][2];
  dist_3b_vec[2][0] = dist_3b[2][0];
  dist_3b_vec[2][1] = dist_3b[2][1];
  dist_3b_vec[2][2] = dist_3b[2][2];
  vector <int> type_3b_vec(3);
  type_3b_vec[0] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype3b[0]));
  type_3b_vec[1] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype3b[1]));
  type_3b_vec[2] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype3b[2]));
  //type_3b_vec[0] = chimes_ptr->atmtoidx[atype3b[0]];
  //type_3b_vec[1] = chimes_ptr->atmtoidx[atype3b[1]];
  //type_3b_vec[2] = chimes_ptr->atmtoidx[atype3b[2]];
  vector<vector<double*> >force_3b_vec;
  force_3b_vec.resize(3,vector<double*>(3));
  force_3b_vec[0][0] = &f3b[0][0];
  force_3b_vec[0][1] = &f3b[0][1];
  force_3b_vec[0][2] = &f3b[0][2];
  force_3b_vec[1][0] = &f3b[1][0];
  force_3b_vec[1][1] = &f3b[1][1];
  force_3b_vec[1][2] = &f3b[1][2];
  force_3b_vec[2][0] = &f3b[2][0];
  force_3b_vec[2][1] = &f3b[2][1];
  force_3b_vec[2][2] = &f3b[2][2];
  vector<double*> stress_vec(9);
  stress_vec[0] = &stress[0];
  stress_vec[1] = &stress[1];
  stress_vec[2] = &stress[2];
  stress_vec[3] = &stress[3];
  stress_vec[4] = &stress[4];
  stress_vec[5] = &stress[5];
  stress_vec[6] = &stress[6];
  stress_vec[7] = &stress[7];
  stress_vec[8] = &stress[8];
  chimes_ptr->compute_3B(dr_3b_vec, dist_3b_vec, type_3b_vec, force_3b_vec, stress_vec, *epot);
  // save forces and stress tensor components 
  f3b[0][0] = *force_3b_vec[0][0];
  f3b[0][1] = *force_3b_vec[0][1];
  f3b[0][2] = *force_3b_vec[0][2];
  f3b[1][0] = *force_3b_vec[1][0];
  f3b[1][1] = *force_3b_vec[1][1];
  f3b[1][2] = *force_3b_vec[1][2];
  f3b[2][0] = *force_3b_vec[2][0];
  f3b[2][1] = *force_3b_vec[2][1];
  f3b[2][2] = *force_3b_vec[2][2];
  stress[0] = *stress_vec[0];
  stress[1] = *stress_vec[1];
  stress[2] = *stress_vec[2];
  stress[3] = *stress_vec[3];
  stress[4] = *stress_vec[4];
  stress[5] = *stress_vec[5];
  stress[6] = *stress_vec[6];
  stress[7] = *stress_vec[7];
  stress[8] = *stress_vec[8];
}

void chimes_compute_4b_props(double dr_4b[6], double dist_4b[6][3], char *atype4b[4], double f4b[4][3], double stress[9], double *epot) {
  // convert all doubles, etc., from C to type vector for C++
  vector <double> dr_4b_vec(6);
  dr_4b_vec[0] = dr_4b[0];
  dr_4b_vec[1] = dr_4b[1];
  dr_4b_vec[2] = dr_4b[2];
  dr_4b_vec[3] = dr_4b[3];
  dr_4b_vec[4] = dr_4b[4];
  dr_4b_vec[5] = dr_4b[5];
  vector< vector<double> > dist_4b_vec(6, vector<double> (3));
  dist_4b_vec[0][0] = dist_4b[0][0];
  dist_4b_vec[0][1] = dist_4b[0][1];
  dist_4b_vec[0][2] = dist_4b[0][2];
  dist_4b_vec[1][0] = dist_4b[1][0];
  dist_4b_vec[1][1] = dist_4b[1][1];
  dist_4b_vec[1][2] = dist_4b[1][2];
  dist_4b_vec[2][0] = dist_4b[2][0];
  dist_4b_vec[2][1] = dist_4b[2][1];
  dist_4b_vec[2][2] = dist_4b[2][2];
  dist_4b_vec[3][0] = dist_4b[3][0];
  dist_4b_vec[3][1] = dist_4b[3][1];
  dist_4b_vec[3][2] = dist_4b[3][2];
  dist_4b_vec[4][0] = dist_4b[4][0];
  dist_4b_vec[4][1] = dist_4b[4][1];
  dist_4b_vec[4][2] = dist_4b[4][2];
  dist_4b_vec[5][0] = dist_4b[5][0];
  dist_4b_vec[5][1] = dist_4b[5][1];
  dist_4b_vec[5][2] = dist_4b[5][2];
  vector <int> type_4b_vec(4);
  type_4b_vec[0] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype4b[0]));
  type_4b_vec[1] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype4b[1]));
  type_4b_vec[2] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype4b[2]));
  type_4b_vec[3] = distance(chimes_ptr->atmtyps.begin(),find(chimes_ptr->atmtyps.begin(), chimes_ptr->atmtyps.end(), atype4b[3]));
  //type_4b_vec[0] = chimes_ptr->atmtoidx[atype4b[0]];
  //type_4b_vec[1] = chimes_ptr->atmtoidx[atype4b[1]];
  //type_4b_vec[2] = chimes_ptr->atmtoidx[atype4b[2]];
  //type_4b_vec[3] = chimes_ptr->atmtoidx[atype4b[3]];
  vector<vector<double*> >force_4b_vec;
  force_4b_vec.resize(4,vector<double*>(3));
  force_4b_vec[0][0] = &f4b[0][0];
  force_4b_vec[0][1] = &f4b[0][1];
  force_4b_vec[0][2] = &f4b[0][2];
  force_4b_vec[1][0] = &f4b[1][0];
  force_4b_vec[1][1] = &f4b[1][1];
  force_4b_vec[1][2] = &f4b[1][2];
  force_4b_vec[2][0] = &f4b[2][0];
  force_4b_vec[2][1] = &f4b[2][1];
  force_4b_vec[2][2] = &f4b[2][2];
  force_4b_vec[3][0] = &f4b[3][0];
  force_4b_vec[3][1] = &f4b[3][1];
  force_4b_vec[3][2] = &f4b[3][2];
  vector<double*> stress_vec(9);
  stress_vec[0] = &stress[0];
  stress_vec[1] = &stress[1];
  stress_vec[2] = &stress[2];
  stress_vec[3] = &stress[3];
  stress_vec[4] = &stress[4];
  stress_vec[5] = &stress[5];
  stress_vec[6] = &stress[6];
  stress_vec[7] = &stress[7];
  stress_vec[8] = &stress[8];
  chimes_ptr->compute_4B(dr_4b_vec, dist_4b_vec, type_4b_vec, force_4b_vec, stress_vec, *epot);
  // save forces and stress tensor components 
  f4b[0][0] = *force_4b_vec[0][0];
  f4b[0][1] = *force_4b_vec[0][1];
  f4b[0][2] = *force_4b_vec[0][2];
  f4b[1][0] = *force_4b_vec[1][0];
  f4b[1][1] = *force_4b_vec[1][1];
  f4b[1][2] = *force_4b_vec[1][2];
  f4b[2][0] = *force_4b_vec[2][0];
  f4b[2][1] = *force_4b_vec[2][1];
  f4b[2][2] = *force_4b_vec[2][2];
  f4b[3][0] = *force_4b_vec[3][0];
  f4b[3][1] = *force_4b_vec[3][1];
  f4b[3][2] = *force_4b_vec[3][2];
  stress[0] = *stress_vec[0];
  stress[1] = *stress_vec[1];
  stress[2] = *stress_vec[2];
  stress[3] = *stress_vec[3];
  stress[4] = *stress_vec[4];
  stress[5] = *stress_vec[5];
  stress[6] = *stress_vec[6];
  stress[7] = *stress_vec[7];
  stress[8] = *stress_vec[8];
}
