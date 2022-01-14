/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<algorithm>
#include<fstream>

using namespace std;

#include "mpi_chimes_interface.h"

// mpi_chimes_interface member functions

void mpi_chimes_interface::init_chimesFF(string chimesFF_paramfile, int my_rank, int my_nprocs)
{
    // Overloaded. Set some MPI stuff and then call the serial version for the rest

    rank   = my_rank;
    nprocs = my_nprocs;

    // Initialize the chimesFF object, read parameters

    init(rank);
    read_parameters(chimesFF_paramfile);
    set_atomtypes(type_list);   
}

mpi_chimes_interface::mpi_chimes_interface(bool small)
{
    
    // For small systems, allow explicit replication prior to ghost atom construction
    // This should ONLY be done for perfectly crystalline systems
    
    allow_replication = small;
    
    // Initialize Pointers, etc for chimes calculator interfacing

    dist_3b.resize(3);
    dist_4b.resize(6);
    
    dr   .resize(3);
    dr_3b.resize(3,std::vector<double>(3));  
    dr_4b.resize(6,std::vector<double>(3));
    
    force_ptr_2b.resize(2,std::vector<double*>(3));
    force_ptr_3b.resize(3,std::vector<double*>(3));
    force_ptr_4b.resize(4,std::vector<double*>(3));
    
    typ_idxs_2b.resize(2);
    typ_idxs_3b.resize(3);
    typ_idxs_4b.resize(4);
    
    max_2b_cut = 0.0;
    max_3b_cut = 0.0;
    max_4b_cut = 0.0;
    
}

mpi_chimes_interface::~mpi_chimes_interface()
{}
    
void mpi_chimes_interface::calculate(int natoms, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps, double & energy, vector<double> & fx, vector<double> & fy, vector<double> & fz, vector<double> & stress)
{    
    // Read system, set up lattice constants/hmats

    // Determine the max outer cutoff (MUST be 2-body, based on ChIMES logic)

    // Initialize private members (LEF)
    max_2b_cut = max_cutoff_2B(true) ;
    max_3b_cut = max_cutoff_3B(true) ;
    max_4b_cut = max_cutoff_4B(true) ;
    
    sys.init(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in, max_2b_cut, allow_replication);    
    
    sys.build_layered_system(atmtyps,poly_orders, max_2b_cut, max_3b_cut, max_4b_cut);

    sys.set_atomtyp_indices(type_list);
    
    sys.run_checks({max_2b_cut,max_3b_cut,max_4b_cut},poly_orders);

    build_neigh_lists(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in);

    
    // Setup vars
    
    int ii, jj, kk, ll, my_start, my_end;
    
    vector<double*> stensor(9);

    for (int idx=0; idx<9; idx++)
        stensor[idx]  = &stress[idx];
    
    ////////////////////////
    // interate over 1- and 2b's 
    ////////////////////////
    
    distribute(my_start, my_end, sys.n_atoms);

    for(int i=my_start; i<=my_end; i++)
    {
        compute_1B(sys.sys_atmtyp_indices[i], energy);

        for(int j=0; j<neighlist_2b[i].size(); j++) // Neighbors of i
        {
            jj = neighlist_2b[i][j];

            dist = sys.get_dist(i,jj,dr); // Populates dr, which is passed by ref (overloaded)
            
            typ_idxs_2b[0] = sys.sys_atmtyp_indices[i ];
            typ_idxs_2b[1] = sys.sys_atmtyp_indices[jj];
        
            force_ptr_2b[0][0] =  &fx[sys.sys_rep_parent[i]];
            force_ptr_2b[0][1] =  &fy[sys.sys_rep_parent[i]];
            force_ptr_2b[0][2] =  &fz[sys.sys_rep_parent[i]];
        
            force_ptr_2b[1][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[jj]]];
            force_ptr_2b[1][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[jj]]];
            force_ptr_2b[1][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[jj]]];        

            compute_2B(dist, dr, typ_idxs_2b, force_ptr_2b, stensor, energy);
        }
    }
    
    ////////////////////////
    // interate over 3b's 
    ////////////////////////

    distribute(my_start, my_end, neighlist_3b.size());
    
    if (poly_orders[1] > 0 )
    {
        for(int i=my_start; i<=my_end; i++)
        {
            ii = neighlist_3b[i][0];
            jj = neighlist_3b[i][1];
            kk = neighlist_3b[i][2];
        
            dist_3b[0] = sys.get_dist(ii,jj,dr_3b[0]); 
            dist_3b[1] = sys.get_dist(ii,kk,dr_3b[1]); 
            dist_3b[2] = sys.get_dist(jj,kk,dr_3b[2]); 
        
            typ_idxs_3b[0] = sys.sys_atmtyp_indices[ii];
            typ_idxs_3b[1] = sys.sys_atmtyp_indices[jj];
            typ_idxs_3b[2] = sys.sys_atmtyp_indices[kk];

            force_ptr_3b[0][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[ii]]];
            force_ptr_3b[0][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[ii]]];
            force_ptr_3b[0][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[ii]]];
        
            force_ptr_3b[1][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[jj]]];
            force_ptr_3b[1][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[jj]]];
            force_ptr_3b[1][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[jj]]];        
        
            force_ptr_3b[2][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[kk]]];
            force_ptr_3b[2][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[kk]]];
            force_ptr_3b[2][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[kk]]];              
        
            compute_3B(dist_3b, dr_3b, typ_idxs_3b, force_ptr_3b, stensor, energy);

        }
    }

    ////////////////////////
    // interate over 4b's 
    ////////////////////////

    distribute(my_start, my_end, neighlist_4b.size());

    if (poly_orders[2] > 0 )
    {
        for(int i=my_start; i<=my_end; i++)
        {
            ii = neighlist_4b[i][0];
            jj = neighlist_4b[i][1];
            kk = neighlist_4b[i][2];
            ll = neighlist_4b[i][3];
        
            dist_4b[0] = sys.get_dist(ii,jj,dr_4b[0]); 
            dist_4b[1] = sys.get_dist(ii,kk,dr_4b[1]); 
            dist_4b[2] = sys.get_dist(ii,ll,dr_4b[2]); 
            dist_4b[3] = sys.get_dist(jj,kk,dr_4b[3]); 
            dist_4b[4] = sys.get_dist(jj,ll,dr_4b[4]); 
            dist_4b[5] = sys.get_dist(kk,ll,dr_4b[5]);         

            typ_idxs_4b[0] = sys.sys_atmtyp_indices[ii];
            typ_idxs_4b[1] = sys.sys_atmtyp_indices[jj];
            typ_idxs_4b[2] = sys.sys_atmtyp_indices[kk];
            typ_idxs_4b[3] = sys.sys_atmtyp_indices[ll];        
        
            force_ptr_4b[0][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[ii]]];
            force_ptr_4b[0][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[ii]]];
            force_ptr_4b[0][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[ii]]];
        
            force_ptr_4b[1][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[jj]]];
            force_ptr_4b[1][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[jj]]];
            force_ptr_4b[1][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[jj]]];        
        
            force_ptr_4b[2][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[kk]]];
            force_ptr_4b[2][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[kk]]];
            force_ptr_4b[2][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[kk]]];      
        
            force_ptr_4b[3][0] =  &fx[sys.sys_rep_parent[sys.sys_parent[ll]]];
            force_ptr_4b[3][1] =  &fy[sys.sys_rep_parent[sys.sys_parent[ll]]];
            force_ptr_4b[3][2] =  &fz[sys.sys_rep_parent[sys.sys_parent[ll]]];                    
              
            compute_4B(dist_4b, dr_4b, typ_idxs_4b, force_ptr_4b, stensor, energy);
        }    
    }

    // Correct for use of replicates, if applicable
    
    energy /= pow(sys.n_replicates+1.0,3.0);
   
    ////////////////////////
    // Finish pressure calculation
    ////////////////////////
    
    for (int idx=0; idx<9; idx++)
        *stensor[idx] /= sys.vol;  
    
}

void mpi_chimes_interface::distribute(int & my_start, int & my_end, int total_items)
{
    int procs_used;

    // Deal with no tasks to perform.
    
    if ( total_items <= 0 ) 
    {
        my_start = 1 ;
        my_end   = 0 ;
        return ;
    }

    // Deal gracefully with more tasks than processors.
    
    if ( nprocs <= total_items ) 
        procs_used = nprocs;
    else
        procs_used = total_items;

    // Use ceil so the last process always has fewer tasks than the other. This improves load balancing.
    
    my_start = ceil( (double) rank * total_items / procs_used);

    if ( rank > total_items ) 
    {
        my_start = total_items + 1;
        my_end = total_items - 1;
    } 
    else if ( rank == procs_used - 1 ) // End of the list.
    {
        my_end = total_items - 1;
    } 
    else // Next starting value - 1 .
    {
        my_end   = ceil( (double) (rank+1) * total_items / procs_used ) - 1;
        if ( my_end > total_items - 1 ) 
            my_end = total_items - 1;
    }
}










