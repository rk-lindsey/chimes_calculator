/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Rebecca K. Lindsey (2020) 
*/

/* ----------------------------------------------------------------------

This code demonstrates how chimesFF{h,cpp} can be used to obtain the 
stress tensor, energy, and per-atom forces for a given system, through the
serial_chimes_interface.

Notes: This script takes as input a standard ChIMES parameter file,
a .xyz file with a, b, and c cell vectors  in the comment line. 
Compile with: 

    g++ -O3 -std=c++11 -o example main.cpp \serial_chimes_interface.cpp \
    chimesFF.cpp
 Run with: 
    ./example <parameter file> <xyz file> 

---------------------------------------------------------------------- */

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>

using namespace std;

#include "serial_chimes_interface.h"    

// Prototypes for some simple helper functions

int    split_line(string line, vector<string> & items);
string get_next_line(istream& str);

int main(int argc, char **argv)
{
    // Read generic code input
    
    string params = argv[1];
    string in_xyz = argv[2];
    
    cout << "Read args:" << " " << params << " " << in_xyz << " " << endl;
    
    // Read the .xyz file    
    
    string            tmp_line;
    vector<string>    tmp_words;
    int                natoms;
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
    if (!coordfile.good())
    {
        cout << "ERROR: Cannot open xyz file " << in_xyz << endl;
        exit(0);
    }
        
    natoms = stoi(get_next_line(coordfile));

    tmp_line = get_next_line(coordfile);
    split_line(tmp_line, tmp_words);

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
        tmp_line = get_next_line(coordfile);
        split_line(tmp_line, tmp_words);
        
        atom_types.push_back(     tmp_words[0] );
        xcrds     .push_back(stod(tmp_words[1]));
        ycrds     .push_back(stod(tmp_words[2]));
        zcrds     .push_back(stod(tmp_words[3]));    
    }
    
    coordfile.close();
    
    // Setup objects to hold the energy, stress tensor, and forces
    
    double                     energy = 0.0;
    vector<double>             stress(9,0.0);   // [xx xy xz yx yy yz zx zy zz]
    vector<vector<double> >    force(natoms);   // [natoms][x, y, or z-component]
    
    for(int i=0; i<natoms; i++)
        force[i].resize(3,0.0);
    
    // Compute ChIMES energy, force, and stress
    
    serial_chimes_interface chimes;        // Create an instance of the serial interface
    
    chimes.init_chimesFF(params, 0);    // Initialize

    chimes.calculate(xcrds, ycrds, zcrds, cell_a, cell_b, cell_c, atom_types, energy, force, stress);


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
        debug_out << scientific << setprecision(6) << force[i][0] << endl 
	              << scientific << setprecision(6) << force[i][1] << endl 
		          << scientific << setprecision(6) << force[i][2] << endl;
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
        cout << "\t" << force[i][0] << "\t" << force[i][1] << "\t" << force[i][2] << endl;

    cout << endl;
}

int split_line(string line, vector<string> & items)
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

string get_next_line(istream& str)
{
    // Read a line and return it, with error checking.
    
    string line;

    getline(str, line);
    
    if ( ! str.good() )
    {
        cout << "Error reading line" << line << endl;
        exit(0);
    }

    return line;
}
