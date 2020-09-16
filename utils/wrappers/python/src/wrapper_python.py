"""

A simple python interface for chimesFF access.

Code Author: Rebecca K. Lindsey (2020)

"""

import ctypes

chimes_wrapper = ctypes.CDLL("libwrapper-C.so")

def get_chimes_max_2b_cutoff():
	return float(chimes_wrapper.get_chimes_max_2b_cutoff())
def get_chimes_max_3b_cutoff():
	return float(chimes_wrapper.get_chimes_max_3b_cutoff())
def get_chimes_max_4b_cutoff():
	return float(chimes_wrapper.get_chimes_max_4b_cutoff())
	
def get_chimes_2b_order():
	return int(chimes_wrapper.get_chimes_2b_order())	
def get_chimes_3b_order():
	return int(chimes_wrapper.get_chimes_3b_order())
def get_chimes_4b_order():
	return int(chimes_wrapper.get_chimes_4b_order())

def set_chimes():
	chimes_wrapper.set_chimes()
	return
	
def init_chimes(rank=0):
	chimes_wrapper.init_chimes(ctypes.int(rank))
	return
	
def read_params(param_file):
	in_param_file = ctypes.c_char_p(param_file.encode())
	chimes_wrapper.chimes_read_params(in_param_file)
	return
	
	
def chimes_compute_2b_props(rij, dr, atype2b, force, stress, epot):
	
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
		ctypes.byref(in_epot)
		)

	return in_force, in_stress, in_epot

def chimes_compute_3b_props(dr_3b, dist_3b, atype3b, force, stress, epot):
	
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
		ctypes.byref(in_epot)
		)

	return in_force, in_stress, in_epot
	

def chimes_compute_4b_props(dr_4b, dist_4b, atype4b, force, stress, epot):
	
	print("in 4b")
	
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
	return in_force, in_stress, in_epot	