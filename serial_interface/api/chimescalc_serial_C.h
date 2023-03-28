/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Nir Goldman (2020) 
*/

#ifdef __cplusplus
extern "C" {
#endif 

void *chimes_open_instance();
void chimes_close_instance(void *handle);
void set_chimes_serial(int small);
void set_chimes_serial_instance(void *handle, int small);
void init_chimes_serial(char *param_file, int *rank);
void init_chimes_serial_instance(void *handle, char *param_file, int *rank);
void calculate_chimes(int natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9]); 
void calculate_chimes_instance(void *handle, int natom, double *xc, double *yc, double *zc, char *atom_types[], double ca[3], double cb[3], double cc[3], double *energy, double fx[], double fy[], double fz[], double stress[9]); 
#ifdef __cplusplus
}
#endif
