/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Rebecca K. Lindsey (2020) 
*/

#ifndef _chimesFF_h
#define _chimesFF_h

#include<vector>
#include<string>
#include<map>

using namespace std;

// Notes:
//
// 1. A Morse-style coordinate transformation is hard-coded (see set_cheby_polys) 
// 2. Polynomials are hard-coded over the domain [-1,1]
// 3. A cubic style cutoff is assumed, and Tersoff is the only other style considered (see get_fcut)

class chimesFF
{
    public:
    
        ////////////////////////
        // General parameters 
        ////////////////////////
	
	int              rank;           // Used to prevent multiple cout statements when accessed from MPI
        int              natmtyps;       // How many atom types are defined for this force field?

        
        vector<int>      poly_orders;    // [bodiedness-1]; i.e. 12 = 2-body only, 12th order; 12 5 = 2+3-body, 0 5 = 3-body only, 5th order
        vector<string>   atmtyps;        // Atom types 
	vector<double>   masses;         // Atom masses
	

        ////////////////////////
        // Functions
        ////////////////////////
    
        chimesFF();
        ~chimesFF();
        
        void init(int mpi_rank);
        
        void read_parameters(string paramfile); 
        
        void compute_1B(const int typ_idx, double & energy );
        
        void compute_2B(const double dx, const vector<double> & dr, const vector<int> typ_idxs, vector<vector<double* > > force, vector<double*> stress, double & energy );
        
        void compute_3B(const vector<double> & dx, const vector<vector<double> > & dr, const vector<int> & typ_idxs, vector<vector<double* > > force, vector<double*> stress, double & energy );

        void compute_4B(const vector<double> & dx, const vector<vector<double> > & dr, const vector<int> & typ_idxs, vector<vector<double* > > force, vector<double*> stress, double & energy );

	void get_cutoff_2B(vector<vector<double> >  & cutoff_2b);	// Populates the 2b cutoffs
	
        double max_cutoff_2B(bool silent = false);    // Returns the largest 2B cutoff
        double max_cutoff_3B(bool silent = false);    // Returns the largest 3B cutoff
        double max_cutoff_4B(bool silent = false);    // Returns the largest 4B cutoff
        
        void set_atomtypes(vector<string> & type_list);
	
	int get_atom_pair_index(int pair_id);

    
    private:
        
        
        string            xform_style;    //  Morse, direct, inverse, etc...
        string            fcut_type;      // cutoff function style (tersoff/cubic)
        double            fcut_var;       // tersoff distance (if fcut_type)
        vector<double>    morse_var;      // [npairs]; morse_lambda
        vector<double>    penalty_params; // [2];  Second dimension: [0] = A_pen, [1] = d_pen
        vector<double>    energy_offsets; // [natmtyps]; Single atom ChIMES energies
        
        // Names (chemical symbols for constituent atoms) .. handled differently for 2-body versus >2-body interactions

        vector<string> pair_params_atm_chem_1;    //[npairs]; // first atom in pair
        vector<string> pair_params_atm_chem_2;    //[npairs]; // second atom in pair
        
        vector<vector<string> > trip_params_atm_chems;    //[ntrips][3]    // Gives chemical symbol  for each ATOM in the triplet (i.e. "Si")    
        vector<vector<string> > trip_params_pair_typs;    //[ntrips][3]    // Gives chemical symbols for each PAIR in the triplet (i.e. "SiO")    
        
        vector<vector<string> > quad_params_atm_chems;    //[quads][3]    // Gives chemical symbol  for each ATOM in the quadruplet (i.e. "Si")    
        vector<vector<string> > quad_params_pair_typs;    //[quads][3]    // Gives chemical symbols for each PAIR in the quadruplet (i.e. "SiO")            

        int n_pair_maps;    // Number of pair maps entries
        int n_trip_maps;    // Number of trip maps entries
        int n_quad_maps;    // Number of quad maps entries
    
        int pair_type_idx;
        int trip_type_idx;
        int quad_type_idx;    
        
        ////////////////////////
        // Definitions for pair, triplet, and quadruplet types
        ////////////////////////

        // 2-body maps

        vector    <string> atom_typ_pair_map; // [nmaps] "slow" maps, based on atom chemical symbol    // Used to build int map -- gives chemical symbol list (i.e. "SiO")
        vector       <int> atom_idx_pair_map; // [nmaps] "slow" maps, based on atom chemical symbol    // Used to build int map -- gives correspoding parameter index (i.e. 5)
        vector       <int> atom_int_pair_map; // [nmaps] "fast" maps, based on atom type index
        vector    <string> atom_int_prpr_map; // [nmaps] "fast" maps, based on atom type index ... returns the "proper" pair type instead of an index
        
        // 3-body maps
        
        vector    <string>    atom_typ_trip_map;    // [nmaps] "slow" maps, based on atom chemical symbol    // Used to build int map -- gives chemical symbol list (i.e. "SiOSiOOO")
        vector       <int>    atom_idx_trip_map;    // [nmaps] "slow" maps, based on atom chemical symbol    // Used to build int map -- gives correspoding parameter index (i.e. 3)
        vector       <int>    atom_int_trip_map;    // [nmaps] "fast" maps, based on atom type index         // gives the correspoding parameter index (i.e. 3) for a unique integer built from type index of three atoms of arbitrary order

        // 4-body maps
        
        vector    <string>    atom_typ_quad_map;    // [nmaps] "slow" maps, based on atom chemical symbol    // Used to build int map -- gives chemical symbol list (i.e. "SiOSiOOO")
        vector      <int>        atom_idx_quad_map; // [nmaps] "slow" maps, based on atom chemical symbol    // Used to build int map -- gives correspoding parameter index (i.e. 3)
        vector      <int>        atom_int_quad_map; // [nmaps] "fast" maps, based on atom type index         // gives the correspoding parameter index (i.e. 3) for a unique integer built from type index of four atoms of arbitrary order
        
        ////////////////////////
        // Polynomial parameters 
        ////////////////////////
    
        // number of coefficients for the pair/triplet/quadruplet type

        vector          <int>   ncoeffs_2b;        // [npairs]

        vector<vector<int>    > chimes_2b_pows;    // [npairs][npowers] power for the coresponding parameter

        vector<vector<double> > chimes_2b_params;    // [npairs][npowers] 2-body polynomial coefficients 
        vector<vector<double> > chimes_2b_cutoff;    // [npairs][2] inner and outer cutoff for pair

        vector<int>                      ncoeffs_3b;          // [ntrips]
        vector<vector<vector<int> > >    chimes_3b_powers;    // [ntrips][nparams][constit. pair]
        vector<vector<double> >          chimes_3b_params;    // [ntrips][nparams]    
        vector<vector<vector<double> > > chimes_3b_cutoff;    // [ntrips][2][constit. pair] inner and outer cutoff for pair 1

        
        vector<int>                      ncoeffs_4b;          // [nquads]
        vector<vector<vector<int> > >    chimes_4b_powers;    // [nquads][nparams][constit. pair]
        vector<vector<double> >          chimes_4b_params;    // [nquads][nparams]    
        vector<vector<vector<double> > > chimes_4b_cutoff;    // [nquads][2][constit. pair] inner and outer cutoff for pair 1

        // Tools for compute functions
        
        inline void set_cheby_polys(double *Tn, double *Tnd, const double dx, const int pair_idx, const double inner_cutoff, const double outer_cutoff, const int bodiedness_idx);
        
        inline void get_fcut(const double dx, const double outer_cutoff, double & fcut, double & fcutderiv);
        
        inline void get_penalty(const double dx, const int & pair_idx, double & E_penalty, double & force_scalar);
        
        inline void build_atom_and_pair_mappers(const int natoms, const int npairs, const vector<int> typ_idxs, const vector<vector<string> > & clu_params_atm_chems, const int & cluidx, vector<int >  & mapped_pair_idx);
        
        int get_proper_pair(string ty1, string ty2);
        
        double max_cutoff(int ntypes, vector<vector<vector<double> > > & cutoff_list);
        
        // Tools for reading the input file
        
        int split_line(string line, vector<string> & items);
        
        string get_next_line(istream& str);        
        
        // Fun stuff
        
        void print_pretty_stuff();
};

#endif
























