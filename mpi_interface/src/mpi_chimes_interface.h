/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

/* ----------------------------------------------------------------------

This class demonstrates how chimesFF{h,cpp} can be used with MPI to 
obtain the stress tensor, energy, and per-atom forces for a given system. 
See main.cpp for a usage example.

Notes: This class has been written for readability rather than speed.
Optimization is recommended prior to use with large systems.

---------------------------------------------------------------------- */

#ifndef _mpi_chimes_interface_h
#define _mpi_chimes_interface_h

#include<vector>
#include<string>

using namespace std;

#include "chimesFF.h"    
#include "serial_chimes_interface.h"


class mpi_chimes_interface : public serial_chimes_interface
{
    public:

        int nprocs = 1;	// Number of processors
        int rank   = 0;	// Index of current processor

        mpi_chimes_interface(bool small = true, int my_rank = 0, int my_nprocs = 1);
        ~mpi_chimes_interface();

        void    calculate(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps, double & energy, vector<double> & fx, vector<double> & fy, vector<double> & fz, vector<double> & stress);
        
    protected:
        
        void distribute(int & my_start, int & my_end, int total_items);
};

#endif






























