"""

    Configuration file for the ChIMES potential energy surface generator (pes_generator.py)

	Don't forget to set PAIR CHEBYSHEV PENALTY SCALING to zero in the parameter file.

"""


CHMS_REPO  = "/Users/lindsey11/My_codes/Forks/chimes_calculator-fork/"

PARAM_FILE = "params.txt"

PAIRTYPES  = [0,    3,    5   ] # Pair type index for scans, i.e. number after "PAIRTYPE PARAMS:" in parameter file
PAIRSTART  = [1.0,  1.0,  1.0 ] # Smallest distance for scan
PAIRSTOP   = [4.0,  4.0,  4.0 ] # Largest distance for scan
PAIRSTEP   = [0.01, 0.01, 0.01] # Step size for scan

TRIPTYPES  = [1,    4   ] # Triplet type index for scans, i.e. number after "TRIPLETTYPE PARAMS:" in parameter file
TRIPSTART  = [1.0,  1.0 ] # Smallest distance for scan
TRIPSTOP   = [4.0,  4.0 ] # Largest distance for scan
TRIPSTEP   = [0.10, 0.10] # Step size for scan

QUADTYPES  = [7   ] # Triplet type index for scans, i.e. number after "TRIPLETTYPE PARAMS:" in parameter file
QUADSTART  = [1.0 ] # Smallest distance for scan
QUADSTOP   = [4.0 ] # Largest distance for scan
QUADSTEP   = [0.10] # Step size for scan
