/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Nir Goldman (2020) 
*/

#ifdef __cplusplus
extern "C" {
#endif 

void set_chimes(int small);
void set_chimes_fromF90(int *small);
void init_chimes(char *param_file, int *rank);
void calculate_chimes(int natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9]); 
void calculate_chimes_fromF90(int *natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9]);
#ifdef __cplusplus
}
#endif

