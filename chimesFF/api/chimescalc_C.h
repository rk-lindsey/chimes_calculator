/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Nir Goldman (2020) 
*/

#ifdef __cplusplus
extern "C" {
#endif 

double get_chimes_max_2b_cutoff();
double get_chimes_max_3b_cutoff();
double get_chimes_max_4b_cutoff();
int get_chimes_2b_order();
int get_chimes_3b_order();
int get_chimes_4b_order();
void set_chimes();
void init_chimes(int rank);
void chimes_read_params(char *param_file);
void chimes_build_pair_int_trip_map();
void chimes_build_pair_int_quad_map();
void chimes_compute_2b_props(double rij, double dr[3], char *atype2b[2], double force[2][3], double stress[9], double *epot);
void chimes_compute_3b_props(double dr_3b[3], double dist_3b[3][3], char *atype3b[3], double f3b[3][3], double stress[9], double *epot);
void chimes_compute_4b_props(double dr_4b[6], double dist_4b[6][3], char *atype4b[4], double f4b[4][3], double stress[9], double *epot);

void chimes_compute_2b_props_fromf90(double *rij, double dr[3], char *type1, char *type2, double force[2][3], double stress[9], double *epot);
void chimes_compute_3b_props_fromf90(double dr_3b[3], double dist_3b[3][3], char *type1, char *type2, char *type3, double f3b[3][3], double stress[9], double *epot);
void chimes_compute_4b_props_fromf90(double dr_4b[6], double dist_4b[6][3], char *type1, char *type2, char *type3, char *type4, double f4b[4][3], double stress[9], double *epot);

#ifdef __cplusplus
}
#endif

