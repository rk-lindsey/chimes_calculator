/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Rebecca K. Lindsey (LLNL)
------------------------------------------------------------------------- */


#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "my_page.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"
#include "pair_chimes.h"
#include "group.h"
#include "update.h"

#include <vector>
#include <iostream>

using namespace LAMMPS_NS;



/*	Functions required by LAMMPS:


settings 	(done)		reads the input script line with arguments defined here
coeff		(done)		set coefficients for one i,j pair type
compute		(done)		workhorse routine that computes pairwise interactions
init_one	(done)		perform initalization for one i,j type pair
init_style 	(done)		initialization specific to this pair style

write_restart			write i,j pair coeffs to restart file
read_restart			read i,j pair coeffs from restart file
write_restart_settings		write global settings to restart file
read_restart_settings		read global settings from restart file
single				force and energy fo a single pairwise interaction between two atoms
*/


PairCHIMES::PairCHIMES(LAMMPS *lmp) : Pair(lmp)
{
	restartinfo = 0;

	int me = comm->me;
	MPI_Comm_rank(world,&me);
	
	chimes_calculator.init(me);
	
	// 2, 3, and 4-body vars for chimesFF access

	dr     .resize(3);
	dr_3b  .resize(3, std::vector<double>(3)); 
	dr_4b  .resize(6, std::vector<double>(3));
	
	dist_3b.resize(3);			       
	dist_4b.resize(6);			       
	
	force_ptr_2b.resize(2,std::vector<double*>(3));
	force_ptr_3b.resize(3,std::vector<double*>(3));
	force_ptr_4b.resize(4,std::vector<double*>(3));
	
	typ_idxs_2b.resize(2);
	typ_idxs_3b.resize(3);
	typ_idxs_4b.resize(4);
	
	// Vars for neighlist construction
	
	tmp_3mer.resize(3);
	tmp_4mer.resize(4);   
	
	if (chimes_calculator.rank == 0)
	{ 
		std::cout << std::endl;
		std::cout << "*********************** WARNING (pair_style chimesFF) ***********************" << std::endl;
		std::cout << "Assuming 2-body interactions have larger cutoffs than (n>2)-body interactions" << std::endl;
		std::cout << "*********************** WARNING (pair_style chimesFF) ***********************" << std::endl;
		std::cout << std::endl;
	}
		
}

PairCHIMES::~PairCHIMES()
{
	if (allocated) 
	{   	    
	    memory->destroy(setflag);
	    memory->destroy(cutsq);
	}
}	

void PairCHIMES::settings(int narg, char **arg)
{
	if (narg != 0) 
		error -> all(FLERR,"Illegal pair_style command. Expects no arguments beyond pair_style name.");

	return;	
}

void PairCHIMES::coeff(int narg, char **arg)
{
	// Expect: pair_coeff * * <parameter file name>

	if (narg != 3) 
		error -> all(FLERR,"Illegal pair_style command. Expects \"pair_coeff * * <parameter file name>\" ");
	
	chimesFF_paramfile = arg[2]; 
	
	chimes_calculator.read_parameters(chimesFF_paramfile);
	
	// Set special LAMMPS flags/cutoffs
	
	if (!allocated)
		allocate();
	
	for(int i=1; i<=chimes_calculator.natmtyps; i++)
	{
		for(int j=i; j<=chimes_calculator.natmtyps; j++)
		{
			setflag[i][j] = 1;
			setflag[j][i] = 1;

			cutsq[i][j]  = chimes_calculator.chimes_2b_cutoff[ chimes_calculator.atom_idx_pair_map[ (i-1)*chimes_calculator.natmtyps + (j-1) ] ][1];
			cutsq[i][j] *= cutsq[i][j];
			
			if (i!=j)
			{
				cutsq[j][i]  = chimes_calculator.chimes_2b_cutoff[ chimes_calculator.atom_idx_pair_map[ (j-1)*chimes_calculator.natmtyps + (i-1) ] ][1];
				cutsq[j][i] *= cutsq[j][i];
			}			
		}
	}
	
	maxcut_3b = chimes_calculator.max_cutoff_3B();
	maxcut_4b = chimes_calculator.max_cutoff_4B();
}

void PairCHIMES::allocate()
{
	allocated = 1;
	
	memory->create(setflag,chimes_calculator.natmtyps+1,chimes_calculator.natmtyps+1,"pair:setflag");
	
	for(int i=1; i<=chimes_calculator.natmtyps; i++)
		for(int j=i; j<=chimes_calculator.natmtyps; j++)
			setflag[i][j] = 0;

	memory->create(cutsq,chimes_calculator.natmtyps+1,chimes_calculator.natmtyps+1,"pair:cutsq");
}	

void PairCHIMES::init_style()
{
	if (atom -> tag_enable == 0)
		error -> all(FLERR,"Pair style ChIMES requires atom IDs");
	
	if (force->newton_pair == 0)
		error->all(FLERR,"Pair style ChIMES requires newton pair on");	
	
	// Set up neighbor lists... borrowing this from pair_airebo:
	// need a full neighbor list, including neighbors of ghosts

	int irequest = neighbor->request(this,instance_me);
	neighbor->requests[irequest]->half = 0;
	neighbor->requests[irequest]->full = 1;
	neighbor->requests[irequest]->ghost = 1;
}


double PairCHIMES::init_one(int i, int j)
{
	// Sets the cutoff for each pair interaction.
	// The maximum of the returned values are used to set outer cutoff for neighbor lists
	// WARNING: This means linking won't work properly if 2-b interactions do not have larger cutoffs than all other
	// higher bodied interactions!!
	
	if (setflag[i][j] == 0) 
		error->all(FLERR,"All pair coeffs are not set");
	
	return sqrt(cutsq[i][j]);
	
}

void PairCHIMES::build_mb_neighlists()
{

	if ( (chimes_calculator.poly_orders[1] == 0) &&  (chimes_calculator.poly_orders[2] == 0))
		return;

	// List gets built based on atoms owned by calling proc. 
	
	neighborlist_3mers.clear();
	neighborlist_4mers.clear();
	
	int i,j,k,l,inum,jnum,knum,lnum, ii, jj, kk, ll;		// Local iterator vars
	int *ilist,*jlist,*klist,*llist, *numneigh,**firstneigh;	// Local neighborlist vars
	tagint 	*tag   = atom -> tag;					// Access to global atom indices (sort of like "parent" indices)
	int     itag, jtag, ktag, ltag;					// holds tags	
	double 	**x    = atom -> x;					// Access to system coordinates
	
	double dist_ij;
	
	////////////////////////////////////////
	// Access to neighbor list vars
	////////////////////////////////////////

	inum       = list -> inum; 		// length of the list
	ilist      = list -> ilist; 		// list of i atoms for which neighbor list exists
	numneigh   = list -> numneigh;		// length of each of the ilist neighbor lists
	firstneigh = list -> firstneigh;	// point to the list of neighbors of i	
	
	for (ii = 0; ii < inum; ii++)	// Loop over real atoms (ai)	
	{

		i     = ilist[ii];		
		itag  = tag[i];			
		jlist = firstneigh[i];		
		jnum  = numneigh[i];	
    	
		for (jj = 0; jj < jnum; jj++)	
		{
			valid_3mer = true;
			valid_4mer = true;		
		
			j     = jlist[jj];	
			jtag  = tag[j];		
			j    &= NEIGHMASK;
			
			if (jtag > itag) 
				continue;
				
			// Check ij distance

			dr[0] = x[j][0] - x[i][0];  
			dr[1] = x[j][1] - x[i][1];
			dr[2] = x[j][2] - x[i][2];

			dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			dist_ij = dist;

			if (dist >= maxcut_3b)
				valid_3mer = false;
				
			if (dist >= maxcut_4b)
				valid_4mer = false;	
				
			if (!valid_3mer && !valid_4mer)
				continue;

			klist = firstneigh[i];		
			knum  = numneigh[i];	
			
			for (kk = 0; kk < knum; kk++)	
			{
				k     = klist[kk];	
				ktag  = tag[k];		
				k    &= NEIGHMASK;
				
				if (j == k)
					continue;
				if (jtag > ktag) 
					continue;
				if (ktag > itag) 
					continue;						
							
 	 			// Check ik distance			

				dr[0] = x[k][0] - x[i][0];  
				dr[1] = x[k][1] - x[i][1];
				dr[2] = x[k][2] - x[i][2];

				dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

				if (dist >= maxcut_3b)
					valid_3mer = false;
					
				if (dist >= maxcut_4b)
					valid_4mer = false;	
					
				if (!valid_3mer && !valid_4mer)
				{
					if(dist_ij < maxcut_3b)
						valid_3mer = true;
					if(dist_ij < maxcut_4b)
						valid_4mer = true;
					continue;
				}		

				// Check jk distance			

				dr[0] = x[k][0] - x[j][0];  
				dr[1] = x[k][1] - x[j][1];
				dr[2] = x[k][2] - x[j][2];

				dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
				
				if (dist > maxcut_3b)
					valid_3mer = false;
					
				if (dist >= maxcut_4b)
					valid_4mer = false;
				
				if (!valid_3mer && !valid_4mer)
				{
					if(dist_ij < maxcut_3b)
						valid_3mer = true;
					if(dist_ij < maxcut_4b)
						valid_4mer = true;
					continue;
				}
				
				// If we're here and valid_3mer == true, then add the triplet to the chimes neigh list        
						
				if (valid_3mer)
				{
					tmp_3mer[0] = i;
					tmp_3mer[1] = j;
					tmp_3mer[2] = k;
					
					neighborlist_3mers.push_back(tmp_3mer);
				}

				if (chimes_calculator.poly_orders[2] == 0)
					continue;
					
				if(!valid_4mer)
				{
					if(dist_ij < maxcut_4b)
						valid_4mer = true;
					continue;
				}
				
				
				llist = firstneigh[i];	
				lnum  = numneigh[i];	
					
				for (ll = 0; ll < lnum; ll++)	
				{
					l     = llist[ll];	
					ltag  = tag[l];		
					l    &= NEIGHMASK;
					
					if ((j == l)||(k == l))
						continue;
					if ((jtag > ltag)||(ktag > ltag)) 
						continue;
					
					if (itag > ltag) 
						continue;
						
					// Check il distance			

					dr[0] = x[l][0] - x[i][0];  
					dr[1] = x[l][1] - x[i][1];
					dr[2] = x[l][2] - x[i][2];

					dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

					if (dist >= maxcut_4b)
						continue;	

					// Check jl distance			
	
					dr[0] = x[l][0] - x[j][0];  
					dr[1] = x[l][1] - x[j][1];
					dr[2] = x[l][2] - x[j][2];
	
					dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
	
					if (dist >= maxcut_4b)

						continue;
								
						
					// Check kl distance			
	
					dr[0] = x[l][0] - x[k][0];  
					dr[1] = x[l][1] - x[k][1];
					dr[2] = x[l][2] - x[k][2];
	
					dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
	
					if (dist >= maxcut_4b)

						continue;
		
	
					// If we're here and valid_4mer == true, then add the quadruplet to the chimes neigh list
							
					if (valid_4mer)
					{
						tmp_4mer[0] = i;
						tmp_4mer[1] = j;
						tmp_4mer[2] = k;
						tmp_4mer[3] = l;
						
						neighborlist_4mers.push_back(tmp_4mer);
					}
				}				
			}
		}
	}
	/*
	std::cout << "3B neighbor list is of length:" << neighborlist_3mers.size() << std::endl;
	std::cout << "4B neighbor list is of length:" << neighborlist_4mers.size() << std::endl;
	*/
}

void PairCHIMES::compute(int eflag, int vflag)
{
	// Vars for access to chimesFF compute_XB functions
	
	std::vector            <double*>  stensor(9);	// pointers to system stress tensor
	std::vector            <double>   stmp(9);	// Dummy var to prevent chimesFF from directly modifying LAMMPS stress
	
	// General LAMMPS compute vars
	
	int 	i,j,k,l,inum,jnum, ii, jj;		// Local iterator vars
	int 	*ilist,*jlist, *numneigh,**firstneigh;	// Local neighborlist vars
	int     idx;

	double 	**x    = atom -> x;		// Access to system coordinates
	double 	**f    = atom -> f;		// Access to system forces
	
	int 	*type  = atom -> type;		// Acces to system atom types (countng starts from 1, chimesFF class expects counting from 0!)
	tagint 	*tag   = atom -> tag;		// Access to global atom indices (sort of like "parent" indices)
	int     itag, jtag, ktag, ltag;		// holds tags
	int 	nlocal = atom -> nlocal;	// Number of real atoms owned by current process .. used used to assure force assignments aren't duplicated
	int 	newton_pair = force -> newton_pair;	// Should f_j be automatically set to -f_i (true) or manually calculated (false)
	double  energy;				// pair energy 

	int me = comm->me;
	MPI_Comm_rank(world,&me);	
	
	// Set up vars controlling if energy/pressure (virial) contributions are computed

	if (eflag || vflag) 
	{
  		ev_setup(eflag,vflag);
	}
	else 
	{
		evflag      = 0;
		vflag_fdotr = 0;
  		vflag_atom  = 0;
	}


	////////////////////////////////////////
	// Access to (2-body) neighbor list vars
	////////////////////////////////////////

	inum       = list -> inum; 		// length of the list
	ilist      = list -> ilist; 		// list of i atoms for which neighbor list exists
	numneigh   = list -> numneigh;		// length of each of the ilist neighbor lists
	firstneigh = list -> firstneigh;	// point to the list of neighbors of i
	
	// Build the ChIMES many-body neighbor lists..
	
	build_mb_neighlists();

//	std::cout << "NEIGHBOR-> AGO: " <<  neighbor->ago << std::endl;


/* Previous approach ... for some reason, this leads to drift in the conserved quantity .... would really be much more efficient if it worked!
	
	// Build the ChIMES many-body neighbor lists.. only do so when LAMMPS neighborlist has been updated
	
	if ( neighbor->ago == 0)
	{
		if (chimes_calculator.rank == 0)
			std::cout << "Updating chimesFF neighbor lists..." << std::endl;
			
		build_mb_neighlists();		
		if (chimes_calculator.rank == 0)
		{
			std::cout << "	Rank 0 3-body list size: " << neighborlist_3mers.size() << std::endl;
			std::cout << "	Rank 0 4-body list size: " << neighborlist_4mers.size() << std::endl;
			std::cout << "	...update complete" << std::endl;
		}
	}
*/	
	
	////////////////////////////////////////
	// Compute 1- and 2-body interactions
	////////////////////////////////////////
	
	for (idx=0; idx<9; idx++)
	{
		stensor[idx]  = &stmp[idx];
		*stensor[idx] = 0;
	}


	for (ii = 0; ii < inum; ii++)				// Loop over the atoms owned by the current process
	{
		i     = ilist[ii];				// Index of the current atom
		itag  = tag[i];					// Get i's global atom index (sort of like its "parent")

		jlist = firstneigh[i];				// Neighborlist for atom i
		jnum  = numneigh[i];				// Number of neighbors of atom i
		
		// First, get the single-atom energy contribution
		
		energy = 0.0;
		
		chimes_calculator.compute_1B(type[i]-1, energy);
		
		if(evflag)
			ev_tally_mb(energy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		
		// Now move on to two-body force, stress, and energy
		
		
		for (jj = 0; jj < jnum; jj++)			// Loop over neighbors of i
		{
			j     = jlist[jj];			// Index of the jj atom
			jtag  = tag[j];				// Get j's global atom index (sort of like its "parent")
			j    &= NEIGHMASK;			// Strip possible extra bits of j
				
			
			if (jtag <= itag) // only allow calculation for j<i, since we've requested a full neighbor list
				continue;
				
			// Get distance using ghost atoms... don't need MIC since we're using ghost atoms
			
			dr[0] = x[j][0] - x[i][0];  
			dr[1] = x[j][1] - x[i][1];
			dr[2] = x[j][2] - x[i][2];
			
			dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			
			typ_idxs_2b[0] = type[i]-1;		// Type (index) of the current atom... subtract 1 to account for chimesFF vs LAMMPS numbering convention
			typ_idxs_2b[1] = type[j]-1;
			
			for (idx=0; idx<3; idx++)
			{
				force_ptr_2b[0][idx] = &f[i][idx];
				force_ptr_2b[1][idx] = &f[j][idx];	
			}
			
			// Do the same for stress tensors
			
			for (idx=0; idx<9; idx++)
				*stensor[idx] = 0;

			energy = 0.0;	
		
			chimes_calculator.compute_2B( dist, dr, typ_idxs_2b, force_ptr_2b, stensor, energy);			

			// "Save"/tally up the energy and stresses to the global virial/energy data objects (see pair.cpp ~ line 1000)
			// Compute pressure, (in contrast to chimes_md) AFTER penalty has been added		
			
			if (evflag)
				ev_tally_mb(energy, *stensor[0], *stensor[1], *stensor[2], *stensor[4], *stensor[5], *stensor[8]);
		}
	}

	if (chimes_calculator.poly_orders[1] > 0)
	{

		////////////////////////////////////////
		// Compute 3-body interactions
		////////////////////////////////////////

		for (ii = 0; ii < neighborlist_3mers.size(); ii++)		
		{
			i     = neighborlist_3mers[ii][0];
			j     = neighborlist_3mers[ii][1];
			k     = neighborlist_3mers[ii][2];

			dr_3b[0][0] = x[j][0] - x[i][0];  
			dr_3b[0][1] = x[j][1] - x[i][1];
			dr_3b[0][2] = x[j][2] - x[i][2];
			
			dr_3b[1][0] = x[k][0] - x[i][0];
			dr_3b[1][1] = x[k][1] - x[i][1];
			dr_3b[1][2] = x[k][2] - x[i][2];
			
			dr_3b[2][0] = x[k][0] - x[j][0];
			dr_3b[2][1] = x[k][1] - x[j][1];
			dr_3b[2][2] = x[k][2] - x[j][2];
			
			dist_3b[0] = sqrt(dr_3b[0][0]*dr_3b[0][0] + dr_3b[0][1]*dr_3b[0][1] + dr_3b[0][2]*dr_3b[0][2]);
			dist_3b[1] = sqrt(dr_3b[1][0]*dr_3b[1][0] + dr_3b[1][1]*dr_3b[1][1] + dr_3b[1][2]*dr_3b[1][2]);
			dist_3b[2] = sqrt(dr_3b[2][0]*dr_3b[2][0] + dr_3b[2][1]*dr_3b[2][1] + dr_3b[2][2]*dr_3b[2][2]);
			
			typ_idxs_3b[0] = type[i]-1;
			typ_idxs_3b[1] = type[j]-1;
			typ_idxs_3b[2] = type[k]-1;

			for (idx=0; idx<3; idx++)
			{
				force_ptr_3b[0][idx] = &f[i][idx];
				force_ptr_3b[1][idx] = &f[j][idx];
				force_ptr_3b[2][idx] = &f[k][idx];
			}

			for (idx=0; idx<9; idx++)
				*stensor[idx] = 0;	
				
			energy = 0;
			
			chimes_calculator.compute_3B( dist_3b, dr_3b, typ_idxs_3b, force_ptr_3b, stensor, energy);				
	
			if (evflag)
				ev_tally_mb(energy, *stensor[0], *stensor[1], *stensor[2], *stensor[4], *stensor[5], *stensor[8]);
		}		
	}

	if (chimes_calculator.poly_orders[2] > 0)
	{
		////////////////////////////////////////
		// Compute 4-body interactions
		////////////////////////////////////////
		
		for (ii = 0; ii < neighborlist_4mers.size(); ii++)		
		{
			i     = neighborlist_4mers[ii][0];
			j     = neighborlist_4mers[ii][1];
			k     = neighborlist_4mers[ii][2];
			l     = neighborlist_4mers[ii][3];			
		
			dr_4b[0][0] = x[j][0] - x[i][0];  
			dr_4b[0][1] = x[j][1] - x[i][1];
			dr_4b[0][2] = x[j][2] - x[i][2];
		
			dr_4b[1][0] = x[k][0] - x[i][0];
			dr_4b[1][1] = x[k][1] - x[i][1];
			dr_4b[1][2] = x[k][2] - x[i][2];
		
			dr_4b[2][0] = x[l][0] - x[i][0];
			dr_4b[2][1] = x[l][1] - x[i][1];
			dr_4b[2][2] = x[l][2] - x[i][2];
		
			dr_4b[3][0] = x[k][0] - x[j][0];
			dr_4b[3][1] = x[k][1] - x[j][1];
			dr_4b[3][2] = x[k][2] - x[j][2];

			dr_4b[4][0] = x[l][0] - x[j][0];
			dr_4b[4][1] = x[l][1] - x[j][1];
			dr_4b[4][2] = x[l][2] - x[j][2];

			dr_4b[5][0] = x[l][0] - x[k][0];
			dr_4b[5][1] = x[l][1] - x[k][1];
			dr_4b[5][2] = x[l][2] - x[k][2];

			for(idx=0;idx<6; idx++)
				dist_4b[idx] = sqrt(dr_4b[idx][0]*dr_4b[idx][0] 
				                  + dr_4b[idx][1]*dr_4b[idx][1] 
					          + dr_4b[idx][2]*dr_4b[idx][2]);					  

			typ_idxs_4b[0] = type[i]-1;
			typ_idxs_4b[1] = type[j]-1;
			typ_idxs_4b[2] = type[k]-1;
			typ_idxs_4b[3] = type[l]-1;

			for (idx=0; idx<3; idx++)
			{
				force_ptr_4b[0][idx] = &f[i][idx];
				force_ptr_4b[1][idx] = &f[j][idx];
				force_ptr_4b[2][idx] = &f[k][idx];
				force_ptr_4b[3][idx] = &f[l][idx];
			}

			for (idx=0; idx<9; idx++)
				*stensor[idx] = 0;	

			energy = 0;	
			
			chimes_calculator.compute_4B( dist_4b, dr_4b, typ_idxs_4b, force_ptr_4b, stensor, energy);

			if (evflag)
				ev_tally_mb(energy, *stensor[0], *stensor[1], *stensor[2], *stensor[4], *stensor[5], *stensor[8]);


		}
	}
	

	
	return;
}
							
void PairCHIMES::write_restart(){}			
void PairCHIMES::read_restart(){}				
void PairCHIMES::write_restart_settings(){}		
void PairCHIMES::read_restart_settings(){}	
void PairCHIMES::single(){}					
