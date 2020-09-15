/* Code author: Nir Goldman (2020) */

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
void chimes_compute_2b_props(double rij, double dr[3], char *atype2b[2], double force[2][3], double stress[9], double *epot);
void chimes_compute_3b_props(double dr_3b[3], double dist_3b[3][3], char *atype3b[3], double f3b[3][3], double stress[9], double *epot);
void chimes_compute_4b_props(double dr_4b[6], double dist_4b[6][3], char *atype4b[4], double f4b[4][3], double stress[9], double *epot);
#ifdef __cplusplus
}
#endif

