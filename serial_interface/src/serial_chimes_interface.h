/* Code author: Rebecca K. Lindsey (2020) */

/* ----------------------------------------------------------------------

This class demonstrates how chimesFF{h,cpp} can be used to obtain the 
stress tensor, energy, and per-atom forces for a given system. See 
main.cpp for a usage example.

Notes: This class has been written for readability rather than speed.
Optimization and parallel distribution is recommended prior to use with 
large systems.

---------------------------------------------------------------------- */

#ifndef _serial_chimes_interface_h
#define _serial_chimes_interface_h

#include<vector>
#include<string>

using namespace std;

#include "chimesFF.h"    


class serial_chimes_interface : public chimesFF
{
    public:
            
        serial_chimes_interface();
        ~serial_chimes_interface();
            
        void    init_chimesFF(string chimesFF_paramfile, int layers, int rank);
        void    calculate(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps, double & energy, vector<vector<double> > & force, vector<double> & stress);

    private:
        
        int n_layers;   // number of replicate layers to make
        int n_atoms;    // number of real atoms
        int n_ghost;    // number of real+ghost atoms
        
        vector<double>    hmat;        // System h-matrix
        vector<double>    invr_hmat;   // Inverse h-matrix
        double            vol;         // System volume
        
        vector<int> sys_parent;        // Index of atom i's (real) parent
        
        vector<string>    type_list;   // A list of possible unique atom types and thier corresponding numerical index, per the parameter file
        
        vector<string>    sys_atmtyps;         // Chemical symbols for all (real+ghost) atoms
        vector<int>        sys_atmtyp_indices; // Atom type indices for all (real+ghost) atoms

        vector<double> sys_x;    // System (i.e. ghost+real) x-coordinates
        vector<double> sys_y;    // System (i.e. ghost+real) y-coordinates
        vector<double> sys_z;    // System (i.e. ghost+real) z-coordinates
        
        double extent_x;    // Length of projection of the rotated cell a onto the x axis
        double extent_y;    // Length of projection of the rotated cell b onto the y axis
        double extent_z;    // Length of projection of the rotated cell c onto the z axis
        
        
        double max_2b_cut;    // Maximum 2-body outer cutoff
        double max_3b_cut;    // Maximum 3-body outer cutoff
        double max_4b_cut;    // Maximum 4-body outer cutoff
        
        vector<vector<int> > neighlist_2b;    // [real atom index][list of real/ghost atom neighbors]
        vector<vector<int> > neighlist_3b;    // [interaction set index][list of 3 atoms within interaction range] -- currently unused
        vector<vector<int> > neighlist_4b;    // [interaction set index][list of 4 atoms within interaction range] -- currently unused    
        
        // Pointers, etc for chimes calculator interfacing (2-body only for now)
        // To set up for many body calculations, see the LAMMPS implementation
        
        double                     dist;
        vector        <double>     dist_3b;
        vector        <double>     dist_4b;
        
        vector        <double>     dr;
        vector<vector<double> >    dr_3b;
        vector<vector<double> >    dr_4b;
        
        vector<vector<double*> >   force_ptr_2b;
        vector<vector<double*> >   force_ptr_3b;
        vector<vector<double*> >   force_ptr_4b;
        
        vector<int>                typ_idxs_2b;
        vector<int>                typ_idxs_3b;
        vector<int>                typ_idxs_4b;
    
        void    reorient_system(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps);
        void    build_layered_system(vector<string> & atmtyps);
        double  get_dist(int i,int j);
        double  get_dist(int i, int j, vector<double> & rij);
        void    build_neigh_lists();
        void    set_atomtyp_indices();    
};

#endif






























