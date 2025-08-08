/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

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

class simulation_system
{
    public: 
        
        simulation_system();
        ~simulation_system();
        
        inline double get_dist(int i,int j, vector<double> & rij);
        inline double get_dist(int i,int j, double* rij);   
        inline double get_dist(int i,int j);
        
        void init(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, double max_2b_cut, bool small = false);
        void set_atomtyp_indices(vector<string> & type_list);
        void copy(simulation_system & to);
        void reorient();
        void build_layered_system(vector<string> & atmtyps, vector<int> & poly_orders, double max_2b_cut, double max_3b_cut, double max_4b_cut);
        void build_neigh_lists(vector<int> & poly_orders, vector<vector<int> > & neighlist_2b, vector<vector<int> > & neighlist_3b, vector<vector<int> > & neighlist_4b, double max_2b_cut, double max_3b_cut, double max_4b_cut);
        void run_checks(const vector<double>& max_cuts, vector<int>&poly_orders);
        
        
        bool allow_replication; // If true, replicates coordinates prior to calculation

        int n_replicates; // number of "real" replicate layers to make
        int n_layers;     // number of ghost layers to make
        int n_atoms;      // number of real atoms
        int n_ghost;      // number of real+ghost atoms
        int n_repl;       // "real" replicate atoms - handling of itty-bitty crystalline systems

		double max_cut;

        vector<int>       sys_atmtyp_indices;   // Atom type indices for all (real+ghost) atoms        
        vector<string>    sys_atmtyps;          // Chemical symbols for all (real+ghost) atoms      
    
        vector<double> sys_x;          // System (i.e. ghost+real) x-coordinates
        vector<double> sys_y;          // System (i.e. ghost+real) y-coordinates
        vector<double> sys_z;          // System (i.e. ghost+real) z-coordinates
        vector<int>    sys_parent;     // Index of atom i's (replicant) parent
        vector<int>    sys_rep_parent; // Replicant's parent atom
       
        double         vol;         // System volume    
        
    private: 
        
        vector<double>    hmat;        // System h-matrix
        vector<double>    invr_hmat;   // Inverse h-matrix
        
        // Used by both sys and neigh coordinates
        
        double latcon_a;
        double latcon_b;
        double latcon_c;
        
        double cell_alpha;
        double cell_beta;
        double cell_gamma;
        
        // only used by neigh coordinates
        
        double extent_x;    // Length of projection of the rotated cell a onto the x axis
        double extent_y;    // Length of projection of the rotated cell b onto the y axis
        double extent_z;    // Length of projection of the rotated cell c onto the z axis   
};


class serial_chimes_interface : public chimesFF
{
    public:
            
        serial_chimes_interface(bool small = true);
        ~serial_chimes_interface();
        
        bool allow_replication; // If true, replicates coordinates prior to calculation
           
        void    init_chimesFF(string chimesFF_paramfile, int rank);
        void    calculate(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps, double & energy, vector<vector<double> > & force, vector<double> & stress);

    private:

        simulation_system sys;      // Input system
        simulation_system neigh;    // Re-oriented ss
        
        vector<string>    type_list;   // A list of possible unique atom types and thier corresponding numerical index, per the parameter file
        
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
        vector        <double>     dr_3b;
        vector        <double>     dr_4b;
        
        vector<vector<double*> >   force_ptr_2b;
        
        vector<int>                typ_idxs_2b;
        vector<int>                typ_idxs_3b;
        vector<int>                typ_idxs_4b;

    
        void build_neigh_lists(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in);
};

inline double simulation_system::get_dist(int i,int j, vector<double> & rij)
{
    /* Orthorhombic way
    rij[0] = sys_x[j] - sys_x[i];
    rij[1] = sys_y[j] - sys_y[i];
    rij[2] = sys_z[j] - sys_z[i];
    */

    // Transform atoms to scaled space

    double inv_ix = invr_hmat[0]*sys_x[i] + invr_hmat[1]*sys_y[i] + invr_hmat[2]*sys_z[i];
    double inv_iy = invr_hmat[3]*sys_x[i] + invr_hmat[4]*sys_y[i] + invr_hmat[5]*sys_z[i];
    double inv_iz = invr_hmat[6]*sys_x[i] + invr_hmat[7]*sys_y[i] + invr_hmat[8]*sys_z[i];

    double inv_jx = invr_hmat[0]*sys_x[j] + invr_hmat[1]*sys_y[j] + invr_hmat[2]*sys_z[j];
    double inv_jy = invr_hmat[3]*sys_x[j] + invr_hmat[4]*sys_y[j] + invr_hmat[5]*sys_z[j];
    double inv_jz = invr_hmat[6]*sys_x[j] + invr_hmat[7]*sys_y[j] + invr_hmat[8]*sys_z[j];

    // Get the distance

    double dx = inv_jx - inv_ix;
    double dy = inv_jy - inv_iy;
    double dz = inv_jz - inv_iz;

    // Convert back to standard units

    rij[0] = hmat[0]*dx + hmat[1]*dy + hmat[2]*dz;
    rij[1] = hmat[3]*dx + hmat[4]*dy + hmat[5]*dz;
    rij[2] = hmat[6]*dx + hmat[7]*dy + hmat[8]*dz;

    return sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);

}


inline double simulation_system::get_dist(int i,int j, double *rij)
{
    /* Orthorhombic way
    rij[0] = sys_x[j] - sys_x[i];
    rij[1] = sys_y[j] - sys_y[i];
    rij[2] = sys_z[j] - sys_z[i];
    */

    // Transform atoms to scaled space

    double inv_ix = invr_hmat[0]*sys_x[i] + invr_hmat[1]*sys_y[i] + invr_hmat[2]*sys_z[i];
    double inv_iy = invr_hmat[3]*sys_x[i] + invr_hmat[4]*sys_y[i] + invr_hmat[5]*sys_z[i];
    double inv_iz = invr_hmat[6]*sys_x[i] + invr_hmat[7]*sys_y[i] + invr_hmat[8]*sys_z[i];

    double inv_jx = invr_hmat[0]*sys_x[j] + invr_hmat[1]*sys_y[j] + invr_hmat[2]*sys_z[j];
    double inv_jy = invr_hmat[3]*sys_x[j] + invr_hmat[4]*sys_y[j] + invr_hmat[5]*sys_z[j];
    double inv_jz = invr_hmat[6]*sys_x[j] + invr_hmat[7]*sys_y[j] + invr_hmat[8]*sys_z[j];

    // Get the distance

    double dx = inv_jx - inv_ix;
    double dy = inv_jy - inv_iy;
    double dz = inv_jz - inv_iz;

    // Convert back to standard units

    rij[0] = hmat[0]*dx + hmat[1]*dy + hmat[2]*dz;
    rij[1] = hmat[3]*dx + hmat[4]*dy + hmat[5]*dz;
    rij[2] = hmat[6]*dx + hmat[7]*dy + hmat[8]*dz;

    return sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);

}

inline double simulation_system::get_dist(int i,int j)
{
    vector<double> rij(3);
    
    return get_dist(i,j,rij);
}

#endif





























