"""

    Helper script to generate potential energy surface scans from ChIMES parameter files
    Expects config.py in the running directory
    Currently supports scan generation for up to four body interactions
    See the sample config.py for more details

    Run with python3.X pes_generator.py

"""



import os
import sys
import numpy as np

sys.path.append(os.path.normpath(os.getcwd()))



# Check that the pes_generator config file is present

if not os.path.exists("config.py"):
    
    print("Error: Cannot find pes_generator config file")
    print("       Searched for: config.py in the current directory")
    exit(0)

# Read the config file

import config



# Check that all necessary chimes_calculator files are present 
 
if not os.path.exists(config.CHMS_REPO + "/chimesFF/api/wrapper_py.py"):
    
    print("Error: Cannot find chimes_calculator direct python wrapper")
    print("       Searched for: " + config.CHMS_REPO + "/chimesFF/api/wrapper_py.py")
    print("       Check CHMS_REPO path provided in config.py. Full paths must be used.")
    exit(0)    
    
if not os.path.exists(config.CHMS_REPO + "/chimesFF/examples/python/lib-C_wrapper-direct_interface.so"):
    
    print("Error: Cannot find chimes_calculator C-wrapper library file")
    print("       Searched for: " + config.CHMS_REPO + "/chimesFF/examples/python/lib-C_wrapper-direct_interface.so")
    print("       Check CHMS_REPO path provided in config.py. Full paths must be used.") 
    print("       Ensure the wrapper file has been compiled (e.g. via ./install from " + config.CHMS_REPO + ")")
    exit(0)

 
 
# Initialize the ChIMES calculator

sys.path.append(config.CHMS_REPO + "/chimesFF/api")

import wrapper_py

wrapper_py.chimes_wrapper = wrapper_py.init_chimes_wrapper(config.CHMS_REPO + "/chimesFF/examples/python/lib-C_wrapper-direct_interface.so")
wrapper_py.set_chimes()
wrapper_py.init_chimes()
wrapper_py.read_params(config.PARAM_FILE)


# Read the atom types 

param_file = open(config.PARAM_FILE)
contents   = param_file.readlines()
param_file.close()

no_atom_types = None
pair_types = []
trip_types = []
quad_types = []

for i in range(len(contents)):
    
    line = contents[i]
    
    if "ATOM TYPES:" in line:
        no_atom_types = int(line.split()[-1])
        print("Expecting ", no_atom_types, "atom types")

    if "PAIRTYPE PARAMS:" in line:
        pair_types.append(line.split()[-2:len(line)])
        print("Found pair type", pair_types[-1], "of index", len(pair_types)-1)

    if "TRIPLETTYPE PARAMS:" in contents[i-1]:
        trip_types.append(line.split()[-3:len(line)])
        print("Found triplet type", trip_types[-1], "of index", len(trip_types)-1)
        
    if "QUADRUPLETYPE PARAMS:" in contents[i-1]:
        quad_types.append(line.split()[-4:len(line)])
        print("Found quad_types type", quad_types[-1], "of index", len(quad_types)-1)        


# Do the 2-body scans

dummy_rij    = [0.0, 0.0, 0.0]
dummy_force  = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
dummy_stress = [0.0]*9

if hasattr(config,'PAIRTYPES'):
    for i in range(len(config.PAIRTYPES)): # Iterate over 2-body types

        
        if len(pair_types)-1 < config.PAIRTYPES[i]:
            print("ERROR: Read",len(pair_types), "2-body types from the parameter file")
            print("       A PES scan was requested for type", int(config.PAIRTYPES[i]))
            print("       Check the config and parameter files")
            exit()        

        print("Writing scanfile for pair type",str(i))

        scanfile = open("chimes_scan_2b.type_" + str(i) + ".dat",'w')
    
        scanfile.write("# " + str(config.PAIRSTART[i]) + " " + str(config.PAIRSTOP[i]) + " " + str(config.PAIRSTEP[i]) + "\n")
    
    
        steps = np.arange(config.PAIRSTART[i], config.PAIRSTOP[i], config.PAIRSTEP[i], dtype=float)
    
        for j in range(steps.size): # r_ij distance
        
            energy = 0.0
        
            dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[j], dummy_rij, pair_types[config.PAIRTYPES[i]], dummy_force, dummy_stress, energy)
        
            scanfile.write(str(steps[j]) + " " + str(energy) + '\n')
        
        scanfile.close()

    



# Do the 3-body scans

dummy_rij    = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
dummy_force  = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
dummy_stress = [0.0]*9

if hasattr(config,'TRIPTYPES'):

    for i in range(len(config.TRIPTYPES)): # Iterate over 3-body types

        if len(trip_types)-1 < config.TRIPTYPES[i]:
            print("ERROR: Read",len(trip_types), "3-body types from the parameter file")
            print("       A PES scan was requested for type", int(config.TRIPTYPES[i]))
            print("       Check the config and parameter files")
            exit()


        print("Writing scanfile for triplet type",str(i))

        scanfile_1 = open("chimes_scan_3b.type_"   + str(i) + ".dat",'w')
        scanfile_2 = open("chimes_scan_2+3b.type_" + str(i) + ".dat",'w')
    
        scanfile_1.write("# " + str(config.TRIPSTART[i]) + " " + str(config.TRIPSTOP[i]) + " " + str(config.TRIPSTEP[i]) + "\n")
        scanfile_2.write("# " + str(config.TRIPSTART[i]) + " " + str(config.TRIPSTOP[i]) + " " + str(config.TRIPSTEP[i]) + "\n")
    
    
        steps = np.arange(config.TRIPSTART[i], config.TRIPSTOP[i], config.TRIPSTEP[i], dtype=float)
    
        for j in range(steps.size):  # r_ij distance
        
            for k in range(steps.size): # r_ik distance
        
                for l in range(steps.size): # r_jk distance
            
                    energy       = 0.0
                    dummy_force  = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] # Need to re-declare because compute_2b changes its dimension
                
                    # Get/write the 3-body only energy
                
                    dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_3b_props([steps[j], steps[k], steps[l]], dummy_rij, trip_types[config.TRIPTYPES[i]], dummy_force, dummy_stress, energy)
        
                    scanfile_1.write(str(steps[j]) + " " + str(steps[k]) + " " + str(steps[l]) + " " + str(energy) + '\n')
                
                    # Add/write the 2-body contributions
        
                    dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[j], [0.0, 0.0, 0.0], [trip_types[config.TRIPTYPES[i]][0], trip_types[config.TRIPTYPES[i]][1]], dummy_force, dummy_stress, energy)
                    dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[k], [0.0, 0.0, 0.0], [trip_types[config.TRIPTYPES[i]][0], trip_types[config.TRIPTYPES[i]][2]], dummy_force, dummy_stress, energy)
                    dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[l], [0.0, 0.0, 0.0], [trip_types[config.TRIPTYPES[i]][1], trip_types[config.TRIPTYPES[i]][2]], dummy_force, dummy_stress, energy)
        
                    scanfile_2.write(str(steps[j]) + " " + str(steps[k]) + " " + str(steps[l]) + " " + str(energy) + '\n')
    
        scanfile_1.close()    
        scanfile_2.close()


# Do the 4-body scans 

dummy_rij    = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
dummy_force  = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
dummy_stress = [0.0]*9

if hasattr(config,'QUADTYPES'):

    for i in range(len(config.QUADTYPES)): # Iterate over 3-body types

        if len(quad_types)-1 < config.QUADTYPES[i]:
            print("ERROR: Read",len(quad_types), "4-body types from the parameter file")
            print("       A PES scan was requested for type", int(config.QUADTYPES[i]))
            print("       Check the config and parameter files")
            exit()

        print("Writing scanfile for quadruplet type",str(i))

        scanfile_1 = open("chimes_scan_4b.type_"     + str(i) + ".dat",'w')
        scanfile_2 = open("chimes_scan_2+3+4b.type_" + str(i) + ".dat",'w')
    
        scanfile_1.write("# " + str(config.QUADSTART[i]) + " " + str(config.QUADSTOP[i]) + " " + str(config.QUADSTEP[i]) + "\n")
        scanfile_2.write("# " + str(config.QUADSTART[i]) + " " + str(config.QUADSTOP[i]) + " " + str(config.QUADSTEP[i]) + "\n")
    
    
        steps = np.arange(config.QUADSTART[i], config.QUADSTOP[i], config.QUADSTEP[i], dtype=float)
    
        for j in range(steps.size):  # r_ij distance
        
            for k in range(steps.size): # r_ik distance
        
                for l in range(steps.size): # r_il distance
            
                    for m in range(steps.size): # r_jk distance
                    
                        for n in range(steps.size): # r_jl distance
                        
                            for o in range(steps.size): # r_kl distance
            
                                energy       = 0.0
                                dummy_force  = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] # Need to re-declare because compute_2b changes its dimension
                
                                # Get/write the 4-body only energy
                
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_4b_props([steps[j], steps[k], steps[l], steps[m], steps[n], steps[o]], 
                                                                        dummy_rij, quad_types[config.QUADTYPES[i]], dummy_force, dummy_stress, energy)
        
                                scanfile_1.write(str(steps[j]) + " " + str(steps[k]) + " " + str(steps[l]) + " " + str(steps[m]) + " " + str(steps[n]) + " " + str(steps[o]) + " " + str(energy) +'\n')
                            
                                # Add/write the 3-body contributions

                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_3b_props([steps[j], steps[k], steps[m]], dummy_rij[0:3], [quad_types[config.QUADTYPES[i]][0],quad_types[config.QUADTYPES[i]][1],quad_types[config.QUADTYPES[i]][2]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_3b_props([steps[j], steps[l], steps[n]], dummy_rij[0:3], [quad_types[config.QUADTYPES[i]][0],quad_types[config.QUADTYPES[i]][1],quad_types[config.QUADTYPES[i]][3]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_3b_props([steps[k], steps[l], steps[n]], dummy_rij[0:3], [quad_types[config.QUADTYPES[i]][0],quad_types[config.QUADTYPES[i]][2],quad_types[config.QUADTYPES[i]][3]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_3b_props([steps[m], steps[o], steps[n]], dummy_rij[0:3], [quad_types[config.QUADTYPES[i]][1],quad_types[config.QUADTYPES[i]][2],quad_types[config.QUADTYPES[i]][3]], dummy_force, dummy_stress, energy)

                                # Add/write the 2-body contributions
        
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[j], [0.0, 0.0, 0.0], [quad_types[config.QUADTYPES[i]][0], quad_types[config.QUADTYPES[i]][1]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[k], [0.0, 0.0, 0.0], [quad_types[config.QUADTYPES[i]][0], quad_types[config.QUADTYPES[i]][2]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[l], [0.0, 0.0, 0.0], [quad_types[config.QUADTYPES[i]][0], quad_types[config.QUADTYPES[i]][3]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[m], [0.0, 0.0, 0.0], [quad_types[config.QUADTYPES[i]][1], quad_types[config.QUADTYPES[i]][2]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[n], [0.0, 0.0, 0.0], [quad_types[config.QUADTYPES[i]][1], quad_types[config.QUADTYPES[i]][3]], dummy_force, dummy_stress, energy)
                                dummy_force, dummy_stress, energy = wrapper_py.chimes_compute_2b_props(steps[o], [0.0, 0.0, 0.0], [quad_types[config.QUADTYPES[i]][2], quad_types[config.QUADTYPES[i]][3]], dummy_force, dummy_stress, energy)
                    
                                scanfile_2.write(str(steps[j]) + " " + str(steps[k]) + " " + str(steps[l]) + " " + str(steps[m]) + " " + str(steps[n]) + " " + str(steps[o]) + " " + str(energy) +'\n')
    
        scanfile_1.close()    
        scanfile_2.close()

    
