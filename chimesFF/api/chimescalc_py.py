"""

	A simple python interface for chimesFF access.

	Expects "libwrapper-C.so" in the same directory as this script


	The following must be included in any python script calling this wrapper:

		import chimescalc_py
		chimescalc_py.chimes_wrapper = chimescalc_py.init_chimes_wrapper("/path/to/lib-C_wrapper-direct_interface.so")
		chimescalc_py.set_chimes()
		chimescalc_py.init_chimes()
		chimescalc_py.read_params("some_parameter_file.txt")

    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author: Rebecca K. Lindsey (2020)

"""

import ctypes


chimes_wrapper = None

def init_chimes_wrapper(lib_name):
	return ctypes.CDLL(lib_name)
def get_chimes_max_2b_cutoff():
	""" Returns maximum 2b cutoff from parameter file """
	return float(chimes_wrapper.get_chimes_max_2b_cutoff())
def get_chimes_max_3b_cutoff():
	""" Returns maximum 3b cutoff from parameter file """
	return float(chimes_wrapper.get_chimes_max_3b_cutoff())
def get_chimes_max_4b_cutoff():
	""" Returns maximum 4b cutoff from parameter file """
	return float(chimes_wrapper.get_chimes_max_4b_cutoff())
	
def get_chimes_2b_order():
	""" Returns 2b polynomial order from parameter file """
	return int(chimes_wrapper.get_chimes_2b_order())	
def get_chimes_3b_order():
	""" Returns 3b polynomial order from parameter file """
	return int(chimes_wrapper.get_chimes_3b_order())
def get_chimes_4b_order():
	""" Returns 4b polynomial order from parameter file """
	return int(chimes_wrapper.get_chimes_4b_order())

def set_chimes():
	""" Instantiates the chimesFF object """
	chimes_wrapper.set_chimes()
	
	# Set data types for various function returns (default is int)
	
	chimes_wrapper.get_chimes_max_2b_cutoff.restype = ctypes.c_double
	chimes_wrapper.get_chimes_max_3b_cutoff.restype = ctypes.c_double
	chimes_wrapper.get_chimes_max_4b_cutoff.restype = ctypes.c_double

	return
	
def init_chimes(rank=0):
	""" 
	Initializes the chimesFF object (sets MPI rank)
	Optionally takes an  MPI rank as input
	"""
	chimes_wrapper.init_chimes(ctypes.c_int(rank))
	return
	
def read_params(param_file):
	""" Reads a ChIMES parameter file """
	in_param_file = ctypes.c_char_p(param_file.encode())
	chimes_wrapper.chimes_read_params(in_param_file)
	return
	
	
def chimes_compute_2b_props(rij, dr, atype2b, force, stress, epot):
	"""
	Compute 2-body contributions to force, stress, and energy for a 
	cluster of two atoms
	
	Inputs:
	rij:     A scalar distance between two atoms i and j
	dr:      The distance vector between i and j (list of length 3)
	atype2b: Atom types for i and j (list of length 2)
	force:   Force compnents for i and j (list of lists w/ lengths [2][3])
	         (Updated by function call)
	stress:  Stress tensor components for system 
	         ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) 
	         (Updated by function call)
	epot:    System potential energy (a scalar; updated by function call)
	
	""" 

	in_rij    =   ctypes.c_double(rij)
	in_dr     =  (ctypes.c_double * 3)(dr[0], dr[1], dr[2])
	in_atype  =  (ctypes.c_char_p * 2)(atype2b[0].encode(), atype2b[1].encode())
	in_force  = ((ctypes.c_double * 3) * 2) ((force[0][0], force[0][1], force[0][2]), 
                                             (force[1][0], force[1][1], force[1][2]))
	in_stress =  (ctypes.c_double * 9)(stress[0], stress[1], stress[2], stress[3], stress[4], stress[5], stress[6], stress[7], stress[8])
	in_epot   =   ctypes.c_double(epot)

	chimes_wrapper.chimes_compute_2b_props(
		in_rij,
		in_dr,
		in_atype,
		in_force,
		in_stress,
		ctypes.byref (in_epot)
		)

	return in_force, in_stress, in_epot.value

def chimes_compute_3b_props(dr_3b, dist_3b, atype3b, force, stress, epot):
	"""
	Compute 3-body contributions to force, stress, and energy for a 
	cluster of three atoms
	
	Inputs:
	dr_3b:   A list of scalar distances between 3 atoms i,j, and k
	         ([r_ij, r_ik, r_jk])
	dist_3b: The distance vectors between i,j, and k
		 ([[dx_ij, dy_ij, dz_ij],[dx_ik, dy_ik, dz_ik], [dx_jk, dy_jk, dz_jk]])
	atype3b: Atom types for i, j and k (list of length 3)
	force:   Force compnents for i,j, and k (list of lists w/ lengths [3][3], i.e. [atom][x/y/z component])
	         (Updated by function call)
	stress:  Stress tensor components for system 
	         ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) 
	         (Updated by function call)
	epot:    System potential energy (a scalar; updated by function call)
	
	""" 	
	in_dr      =  (ctypes.c_double * 3)(dr_3b[0], dr_3b[1], dr_3b[2])
	in_dist    = ((ctypes.c_double * 3) * 3)((dist_3b[0][0], dist_3b[0][1], dist_3b[0][2]), 
	                                         (dist_3b[1][0], dist_3b[1][1], dist_3b[1][2]), 
						 (dist_3b[2][0], dist_3b[2][1], dist_3b[2][2]))
	in_atype   =  (ctypes.c_char_p * 3)(atype3b[0].encode(), atype3b[1].encode(), atype3b[2].encode())
	in_force   = ((ctypes.c_double * 3) * 3) ((force[0][0], force[0][1], force[0][2]), 
						  (force[1][0], force[1][1], force[1][2]),
						  (force[2][0], force[2][1], force[2][2]))
	in_stress  =  (ctypes.c_double * 9)(stress[0], stress[1], stress[2], stress[3], stress[4], stress[5], stress[6], stress[7], stress[8])
	in_epot    =   ctypes.c_double(epot)
	
	chimes_wrapper.chimes_compute_3b_props(
		in_dr,
		in_dist,
		in_atype,
		in_force,
		in_stress,
		ctypes.byref (in_epot)
		)

	return in_force, in_stress, in_epot.value
	

def chimes_compute_4b_props(dr_4b, dist_4b, atype4b, force, stress, epot):
	"""
	Compute 4-body contributions to force, stress, and energy for a 
	cluster of four atoms
	
	Inputs:
	dr_4b:   A list of scalar distances between 4 atoms i,j,k, and l
	         ([r_ij, r_ik, r_il, r_jk, r_jl, r_kl])
	dist_4b: The distance vectors between i,j, and k
		 ([[dx_ij, dy_ij, dz_ij],[dx_ik, dy_ik, dz_ik], [dx_il, dy_il, dz_il]],...)
	atype4b: Atom types for i, j, k, and l (list of length 4)
	force:   Force compnents for i,j, k, and l (list of lists w/ lengths [4][3], i.e. [atom][x/y/z component])
	         (Updated by function call)
	stress:  Stress tensor components for sys3tem 
	         ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) 
	         (Updated by function call)
	epot:    System potential energy (a scalar; updated by function call)
	
	""" 	
	in_dr      =  (ctypes.c_double * 6)(dr_4b[0], dr_4b[1], dr_4b[2], dr_4b[3], dr_4b[4], dr_4b[5])
	in_dist    = ((ctypes.c_double * 3) * 6) ((dist_4b[0][0], dist_4b[0][1], dist_4b[0][2]), 
                                              (dist_4b[1][0], dist_4b[1][1], dist_4b[1][2]), 
                                              (dist_4b[2][0], dist_4b[2][1], dist_4b[2][2]), 
                                              (dist_4b[3][0], dist_4b[3][1], dist_4b[3][2]), 
                                              (dist_4b[4][0], dist_4b[4][1], dist_4b[4][2]), 
                                              (dist_4b[5][0], dist_4b[5][1], dist_4b[5][2]))
	in_atype   =  (ctypes.c_char_p * 4)(atype4b[0].encode(), atype4b[1].encode(), atype4b[2].encode(), atype4b[3].encode())
	in_force   = ((ctypes.c_double * 3) * 4) ((force[0][0], force[0][1], force[0][2]), 
                                              (force[1][0], force[1][1], force[1][2]),
                                              (force[2][0], force[2][1], force[2][2]),
                                              (force[3][0], force[3][1], force[3][2]))
	in_stress  =  (ctypes.c_double * 9)(stress[0], stress[1], stress[2], stress[3], stress[4], stress[5], stress[6], stress[7], stress[8])
	in_epot    =   ctypes.c_double(epot)
	
	chimes_wrapper.chimes_compute_4b_props(
		in_dr,
		in_dist,
		in_atype,
		in_force,
		in_stress,
		ctypes.byref(in_epot)
		)
	return in_force, in_stress, in_epot.value	
