/*
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Rebecca K. Lindsey (2022)
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

#include "mpi_chimes_interface.h"
#include "chimescalc_mpi_C.h"

static mpi_chimes_interface chimes; 
static mpi_chimes_interface *chimes_ptr;

void set_chimes_mpi(int small=1)
{
	if ((small!=0)&&(small!=1))
	{
		cout << "ERROR: Small must be set to 0 (false) or 1 (true)" << endl;
		cout << "Received: " << small << endl;
		exit(0);
	}
    
	chimes_ptr = &chimes;
	chimes_ptr->allow_replication = small;
}

void init_chimes_mpi(char *param_file, int *rank, int *nprocs)
{
    chimes_ptr->init_chimesFF(param_file, *rank, *nprocs);
}

void calculate_chimes_mpi(int natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9])
{
    vector<double>    x_vec(natom);
    vector<double>    y_vec(natom);
    vector<double>    z_vec(natom);

    vector<double> fx_vec(natom);
    vector<double> fy_vec(natom);
    vector<double> fz_vec(natom);
    
    vector<string> atom_types_vec;
    atom_types_vec.resize(natom);

    for (int i = 0; i < natom; i++) 
    {
        x_vec[i] = xc[i];
        y_vec[i] = yc[i];
        z_vec[i] = zc[i];
        
        fx_vec[i] = fx[i];
        fy_vec[i] = fy[i];
        fz_vec[i] = fz[i];
        
        atom_types_vec[i] = atom_types[i];
    }
    vector<double> stress_vec(9,0.0);

    for (int i = 0; i < 9; i++) 
        stress_vec[i] = stress[i];

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
    
    chimes_ptr->calculate(natom, x_vec, y_vec, z_vec, cell_a_vec, cell_b_vec, cell_c_vec, atom_types_vec, *energy, fx_vec, fy_vec, fz_vec, stress_vec);
    
    for (int i = 0; i < natom; i++) 
    {
      fx[i] = fx_vec[i];
      fy[i] = fy_vec[i];
      fz[i] = fz_vec[i];
    }
    
    for (int i = 0; i < 9; i++)
        stress[i] = stress_vec[i];

}



