#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define square(x) ((x)*(x))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define MaxPairsPerAtom (long) 90 // maximal number of neighbors in Verlet list

extern long Nall, npm, lpm, lpm1, nrun;
extern double *x, *y, *z, *vx, *vy, *vz, *xold, *yold, *zold; 
extern double *cons_fx, *cons_fy, *cons_fz, *fx_aux, *fy_aux, *fz_aux;
extern double *gauss, *f_aux, *f_cons;

extern long Morse, WCA, LJ, KG;
extern long seed, steps, step_meas, Numb_meas, t_meas;

extern long *List, *Advance, *Marker1, *Marker2;
extern long ListUpdateRequested, MaxListLength;
extern double *DisplaceList_x, *DisplaceList_y, *DisplaceList_z;
extern double skin, Range_Verlet;

extern double e_LJ, LJ_shift, f_ext;
extern double e_M, k_spring, R0, alpha, r_min; // Morse parameters
extern double sigma, Range_WCA, Range_WCA2, Range_LJ, Range_LJ2; 
extern double dt, T, mass, Gamma, inv_mass, inv_sigma2, sqr_Gamma_dt;
extern float ran1();
extern float gasdev();

//New

extern double sigma_b, e_b, Range_LJ_b, inv_sigma_b2;
extern long gen_init, time_width;
extern double count, count2, E, *index_ads, *index_ads_energy, Rg_x, Rg_y, Rg_z, count_sig, count_en;
extern double *Zcor, *Zcor2, *Mcor, *Mcor2;
extern float t_meas_temp, z_end_temp, Rg_x_temp, Rg_y_temp, Rg_z_temp, count_sig_temp, count_en_temp;
extern double *corrChain_z, *corrChain_nrg;
extern int buffer;
//

extern char   fcfg[60], fin1[60], fin2[60], fout1[60], fout2[60];

long  get_params(char *fname);
long  get_monomers(char *);
long  put_monomers(char *);
long  init_mem(void);

void velocity_Verlet(void);
void Verlet_List(void);
void cons_force(void);
void calculate_kinetic_energy(void);
void calculate_potential_energy(long);
void f_list(void);

long put_measurment(char *);
void init_measure(void);
void initialize(void);
void step_measure(long);
void finish(void);


double force_M(double);
double force_KG(double);
double force_LJ(double);
double force_WCA(double);
double Morse_pot(double);
double KG_pot(double);
double LJ_pot(double);
double WCA_pot(double);
void get_gasdev(long);
void check_bonds(void);

//my additions
extern FILE *anim;
extern int nrg_step;
extern double nrg_ncrmnt;
extern long need;


//New
double force_LJ_boundary(double);
double LJ_pot_boundary(double);
void adsorption_measurement(void);
void gyration_radius(void);
void Corr_Along_Chain(void);
//

