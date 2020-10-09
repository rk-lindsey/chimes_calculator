"""

	A simple python interface for the serial_chimes_interface.

	Expects "libwrapper-C.so" in the same directory as this script


	The following must be included in any python script calling this wrapper:

		import wrapper_py
		wrapper_py.chimes_wrapper = wrapper_py.init_chimes_wrapper("/path/to/libwrapper-C.so")
		wrapper_py.set_chimes()
		wrapper_py.init_chimes()
		wrapper_py.read_params("some_parameter_file.txt", <nlayers>)

	Code Author: Rebecca K. Lindsey (2020)

"""

import ctypes


chimes_wrapper = None

def init_chimes_wrapper(lib_name):
	return ctypes.CDLL(lib_name)

def set_chimes():
	""" Instantiates the chimesFF object """
	chimes_wrapper.set_chimes()
	return


def init_chimes(param_file, nlayer, rank):
	""" 
	Initializes the chimesFF object (sets MPI rank)
	Optionally takes an  MPI rank as input
	"""
	in_paramfile = ctypes.c_char_p(param_file.encode())
	in_nlayer    = ctypes.c_int(nlayer)
	in_rank      = ctypes.c_int(rank)
	chimes_wrapper.init_chimes(in_paramfile, ctypes.byref(in_nlayer), ctypes.byref(in_rank))
	return

def calculate_chimes(natoms,xcrd,ycrd,zcrd,atmtyps,cell_a,cell_b,cell_c,energy,fx,fy,fz,stress):
	""" 
	Computes the ChIMES forces, energy, and stress tensor for a given system
	
	Inputs:	
	natoms:	Number of atoms in system 
	xcrd:	System x-coordinates
	ycrd:	System y-coordinates
	zcrd:	System z-coordinates
	atmtyps:System atom types
	cell_a:	System a lattice vector
	cell_b:	System b lattice vector
	cell_c:	System c lattice vector
	energy:	System energy
	fx:	X force components for system atoms
	fy:	Y force components for system atoms
	fz:	Z force components for system atoms
	stress:	System stress tensor
	
	Returns updated fx, fy, fz, stress, and energy
	
	"""
	
	encoded = []
	
	for i in range(natoms):
		encoded.append(atmtyps[i].encode())

	in_natom   = ctypes.c_int(natoms) 
	in_xcrd    = (ctypes.c_double * natoms) (*xcrd)
	in_ycrd    = (ctypes.c_double * natoms) (*ycrd)
	in_zcrd    = (ctypes.c_double * natoms) (*zcrd)
	in_atmtyps = (ctypes.c_char_p * natoms) (*encoded)
	in_cell_a  = (ctypes.c_double * 3) (*cell_a)
	in_cell_b  = (ctypes.c_double * 3) (*cell_b)
	in_cell_c  = (ctypes.c_double * 3) (*cell_c)
	in_energy  = ctypes.c_double(energy)
	in_fx      = (ctypes.c_double * natoms) (*fx)
	in_fy      = (ctypes.c_double * natoms) (*fy)
	in_fz      = (ctypes.c_double * natoms) (*fz)
	in_stress  = (ctypes.c_double * 9) (*stress)

	chimes_wrapper.calculate_chimes(in_natom,
					in_xcrd,   
					in_ycrd,   
					in_zcrd,   
					in_atmtyps,   
					in_cell_a,
					in_cell_b, 
					in_cell_c, 
					ctypes.byref (in_energy), 
					in_fx,     
					in_fy,     
					in_fz,     
					in_stress)
					
	return in_fx, in_fy, in_fz, in_stress, in_energy.value


