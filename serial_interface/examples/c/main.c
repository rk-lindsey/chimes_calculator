/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Nir Goldman (2020) 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "chimescalc_serial_C.h"

int main (int argc, char **argv) 
{

  
  set_chimes(small);
  
    init_chimes(argv[1],&rank);
  
  calculate_chimes(natom, xc, yc, zc, atom, ca, cb, cc, &energy, fx, fy, fz, stress);
  for (i = 0; i < 9; i++) {
    stress[i] *= GPa;
  }
  
  
}
