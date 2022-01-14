/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

/* ----------------------------------------------------------------------

This code demonstrates how chimesFF{h,cpp} can be used to obtain the 
stress tensor, energy, and per-atom forces for a given system, through the
mpi_chimes_interface.

Notes: This script takes as input a standard ChIMES parameter file,
a .xyz file with a, b, and c cell vectors  in the comment line. 
Compile with: 

 See Makefile for compilation instructions.
 
 Run with: 
    srun -n <nproc> ./example <parameter file> <xyz file> <allow replictes (0/1 or true/false)>

---------------------------------------------------------------------- */

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cstring>

using namespace std;

#include "mpi_chimes_interface.h"    

// Setup for use with MPI

#ifdef USE_MPI
    #include <mpi.h>
#endif



// Prototypes for some simple helper functions

int    split_line(string line, vector<string> & items, int rank);
string get_next_line(istream& str, int rank);
void   tally_FES(double & energy,  vector<double> & fx, vector<double> & fy, vector<double> & fz, vector<double> & stress, int rank);

int main(int argc, char **argv)
{
    // Prepare for use with MPI

    int nprocs, rank;

    #ifdef USE_MPI
        MPI_Init     (&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
        if (rank==0)
            cout << "Code compiled in MPI mode."; 
    #else
        if (rank==0)
            cout << "Code compiled in serial mode."; 
    #endif

    if (rank==0)
        cout <<" Will run on " << nprocs << " processor(s)." << endl;

    // Read generic code input
    
    string params = argv[1];
    string in_xyz = argv[2];
    
    bool   is_small = false;
    
    if (rank == 0)
    {
        cout << "Read args:" << endl;
    
        for (int i=1; i<argc; i++)
    
            cout << i << " " << argv[i] << endl;
    }
    
    if(argc == 4)
        if((strncmp(argv[3],"true",4) == 0) || (strncmp(argv[3],"True",4) == 0) || (strncmp(argv[3],"TRUE",4) == 0) || (strncmp(argv[3],"1"   ,1) == 0))
            is_small = true;
    
    
    // Read the .xyz file    
    
    string            tmp_line;
    vector<string>    tmp_words;
    int               natoms;
    double            lx, ly, lz;
    
    vector<double>cell_a(3);
    vector<double>cell_b(3);
    vector<double>cell_c(3);
    
    vector<string>    atom_types;
    vector<double>    xcrds;
    vector<double>    ycrds;
    vector<double>    zcrds;
    
    ifstream coordfile;
    coordfile.open(in_xyz);
    if ((!coordfile.good()) && (rank == 0))
    {
        cout << "ERROR: Cannot open xyz file " << in_xyz << endl;
        exit(0);
    }
        
    natoms = stoi(get_next_line(coordfile, rank));

    tmp_line = get_next_line(coordfile, rank);
    split_line(tmp_line, tmp_words, rank);

    cell_a[0] = stod(tmp_words[0]);
    cell_a[1] = stod(tmp_words[1]);
    cell_a[2] = stod(tmp_words[2]);

    cell_b[0] = stod(tmp_words[3]);
    cell_b[1] = stod(tmp_words[4]);
    cell_b[2] = stod(tmp_words[5]);
    
    cell_c[0] = stod(tmp_words[6]);
    cell_c[1] = stod(tmp_words[7]);
    cell_c[2] = stod(tmp_words[8]);    

    for(int i=0; i<natoms; i++)
    {
        tmp_line = get_next_line(coordfile, rank);
        split_line(tmp_line, tmp_words, rank);
        
        atom_types.push_back(     tmp_words[0] );
        xcrds     .push_back(stod(tmp_words[1]));
        ycrds     .push_back(stod(tmp_words[2]));
        zcrds     .push_back(stod(tmp_words[3]));    
    }
    
    coordfile.close();
    
    // Setup objects to hold the energy, stress tensor, and forces
    
    double                     energy = 0.0;
    vector<double>             stress(9,0.0);   // [xx xy xz yx yy yz zx zy zz]
    vector<double>             fx(natoms,0.0);
    vector<double>             fy(natoms,0.0);
    vector<double>             fz(natoms,0.0);        
    
    // Compute ChIMES energy, force, and stress
    
    mpi_chimes_interface chimes(is_small);        // Create an instance of the MPI interface

    chimes.init_chimesFF(params, rank, nprocs);    // Initialize

    chimes.calculate(natoms, xcrds, ycrds, zcrds, cell_a, cell_b, cell_c, atom_types, energy, fx, fy, fz, stress);
    
    // Now MPI reduce
    
    tally_FES(energy, fx, fy, fz, stress, rank);    

    if(rank == 0)
    {
        #if DEBUG==1
    
            ofstream debug_out;
            debug_out.open("debug.dat");
   
            debug_out << fixed << setprecision(6) << energy << endl;
            debug_out << fixed << setprecision(6) << stress[0]*6.9479 << endl;
            debug_out << fixed << setprecision(6) << stress[4]*6.9479 << endl;
            debug_out << fixed << setprecision(6) << stress[8]*6.9479 << endl;
            debug_out << fixed << setprecision(6) << stress[1]*6.9479 << endl;
            debug_out << fixed << setprecision(6) << stress[2]*6.9479 << endl;
            debug_out << fixed << setprecision(6) << stress[5]*6.9479 << endl;

            for(int i=0; i<natoms; i++)
                debug_out 
              << scientific << setprecision(6) << fx[i] << endl 
                  << scientific << setprecision(6) << fy[i] << endl 
              << scientific << setprecision(6) << fz[i] << endl;
        #endif

        cout << endl;
        cout << "Success! " << endl;
        cout << "Energy (kcal/mol):    " << endl << "\t" << energy << endl;

        cout << "Stress tensors (GPa): " << endl;
        cout << "\ts_xx: " << stress[0]*6.9479 << endl;
        cout << "\ts_yy: " << stress[4]*6.9479 << endl;
        cout << "\ts_zz: " << stress[8]*6.9479 << endl;
        cout << "\ts_xy: " << stress[1]*6.9479 << endl;
        cout << "\ts_xz: " << stress[2]*6.9479 << endl;
        cout << "\ts_yz: " << stress[5]*6.9479 << endl;
        
        cout << "Forces (kcal/mol/A): " << endl;
        for(int i=0; i<natoms; i++)
            cout << "\t" << fx[i] << "\t" << fy[i] << "\t" << fz[i] << endl;

        cout << endl;
    }

#ifdef USE_MPI
	MPI_Finalize() ;
#endif
	exit(0);
}

int split_line(string line, vector<string> & items, int rank)
{
    // Break a line up into tokens based on space separators.
    // Returns the number of tokens parsed.
    
    string       contents;
    stringstream sstream;

    // Strip comments beginining with ! or ## and terminal new line

    int pos = line.find('!');
      
    if ( pos != string::npos ) 
        line.erase(pos, line.length() - pos);

    pos = line.find("##");
    if ( pos != string::npos ) 
        line.erase(pos, line.length()-pos);

    pos = line.find('\n');
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    sstream.str(line);
     
    items.clear();

    while ( sstream >> contents ) 
        items.push_back(contents);

    return items.size();
}

string get_next_line(istream& str, int rank)
{
    // Read a line and return it, with error checking.
    
    string line;

    getline(str, line);
    
    if ( (!str.good()) && (rank == 0))
    {
        cout << "Error reading line" << line << endl;
        exit(0);
    }

    return line;
}

void tally_FES(double & energy,  vector<double> & fx,   vector<double> & fy, vector<double> & fz, vector<double> & stress, int rank)
{
    #ifndef USE_MPI
        return;
    #endif

    // Flatten out force and stress arrays
    
    double *flattened_fx      = (double *) fx.data();
    double *flattened_fy      = (double *) fy.data();
    double *flattened_fz      = (double *) fz.data();
    double *flattened_stress  = (double *) stress.data();    
    int natoms = fx.size();

    // Tally up
    
    MPI_Allreduce(MPI_IN_PLACE, &energy,          1,      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, flattened_fx,     natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, flattened_fy,     natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, flattened_fz,     natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, flattened_stress, 9,      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
}
