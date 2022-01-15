/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Rebecca K. Lindsey (2022) 

    Notes: chimescalc_mpi inherits chimesFF but overloads calculate(...) and init_chimes_serial(...)
*/

#ifdef __cplusplus
extern "C" {
#endif 

void set_chimes_mpi(int small);
void init_chimes_mpi(char *param_file, int *rank, int *nprocs);
void calculate_chimes_mpi(int natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9]); 

#ifdef __cplusplus
}
#endif
