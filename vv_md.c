#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "defs.h"


#define CLOCK_PER_SEC (double) 1000000

long Nall, npm, lpm, lpm1, nrun;
double *x, *y, *z, *vx, *vy, *vz, *xold, *yold, *zold;
double *cons_fx, *cons_fy, *cons_fz, *fx_aux, *fy_aux, *fz_aux;
double *gauss, *f_cons, *f_aux, f_ext;

long seed, Morse, WCA, LJ, KG, steps, step_meas, Numb_meas, t_meas;
long *List, *Advance, *Marker1, *Marker2, ListUpdateRequested, MaxListLength;
double *DisplaceList_x, *DisplaceList_y, *DisplaceList_z, skin;
double sigma, e_LJ, LJ_shift, k_spring, R0, alpha, e_M, r_min;
double dt, T, mass, inv_mass, inv_sigma2, Gamma, sqr_Gamma_dt;
double Range_Verlet, Range_WCA, Range_WCA2, Range_LJ, Range_LJ2;

//New

double sigma_b, e_b, Range_LJ_b, inv_sigma_b2;
long gen_init, time_width;
double count, count2, E, *index_ads, *index_ads_energy, Rg_x, Rg_y, Rg_z, count_sig, count_en;
double *Zcor, *Zcor2, *Mcor, *Mcor2;
float t_meas_temp, z_end_temp, Rg_x_temp, Rg_y_temp, Rg_z_temp, count_sig_temp, count_en_temp;
double *corrChain_z, *corrChain_nrg;
int buffer;

//

float ran1();
char   fcfg[60], fin1[60], fin2[60], fout1[60], fout2[60];

// my additions
FILE *anim;
int nrg_step;
double nrg_ncrmnt;
long need;

int main(long argc, char *argv[])
{	
	
	char *ptr;
	long simID;
	need = strtol(argv[1], &ptr, 10);
	simID = strtol(argv[2], &ptr, 10);
	
	
	FILE *meas_bin_file;
	FILE *meas_bin_file_2;
	char file1[60];
	char file2[60];
	
	for(nrg_step=0;nrg_step<1;nrg_step++){
		clock_t t0, t1;
		long i, j, cnt;
		double dtime, random, noise, half=0.5; 
		FILE *fp, *gp;

		strcpy(fcfg, "in");

		initialize();
		init_measure();

		// sprintf(file1,"Measurments_bin_%lg_%ld",e_b,simID);
		// meas_bin_file = fopen(file1,"wb");
		
		sprintf(file1,"Measurments_bin_Z_%lg_%ld_%ld",e_b,simID,Nall);
		meas_bin_file = fopen(file1,"wb");
		
		sprintf(file2,"Measurments_bin_NRG_%lg_%ld_%ld",e_b,simID,Nall);
		meas_bin_file_2 = fopen(file2,"wb");


		printf(" Let's start...\n");
		t0 = clock();

		// perform initial measurement
		
		// perform the very first simulation step!!!
		Verlet_List(); // 1_st update Verlet list 
		ListUpdateRequested = 0;

		for(j=0;j<Nall;j++) { // remember current monomer coordinates
			xold[j] = x[j]; yold[j] = y[j]; zold[j] = z[j];
		}
		//****************************************************

		cons_force(); // compute initial forces

		for(j=0;j<Nall;j++) {

			random = ran1(&seed);
			noise = sqr_Gamma_dt * (random-half);
			fx_aux[j] = cons_fx[j] - mass*Gamma*vx[j] + noise;

			random = ran1(&seed);
			noise = sqr_Gamma_dt * (random-half);
			fy_aux[j] = cons_fy[j] - mass*Gamma*vy[j] + noise;

			random = ran1(&seed);
			noise = sqr_Gamma_dt * (random-half);
			fz_aux[j] = cons_fz[j] - mass*Gamma*vz[j] + noise;
		}
		//***************************************************
		// start the simulation
		//***************************************************
		
		// loop over total integration steps x runs
  
		//anim = fopen("anim","a");
		
		
		int which_meas = 1;
		sprintf(fout1,"initial_conditions/P_end_%lg_%ld",e_b,Nall);
		
		int k;
		for(k=0;k<=(Nall-buffer)*Numb_meas;k++){
			corrChain_z[k] = 0.0;
			corrChain_nrg[k] = 0.0;
		}
		
		for (t_meas=1;t_meas<=(steps*nrun);t_meas++) {

			velocity_Verlet();
			//check_bonds();
			
			if (!(t_meas%step_meas)){

				Corr_Along_Chain();

				//fprintf(meas_file,"%d,%lg,%lg,%lg\n",t_meas,count,Rg_par,z[Nall-1]);
				//fprintf(meas_txt_file,"%d,%lg,%lg,%lg,%lg,%lg,%lg\n",t_meas,count_sig,count_en,Rg_x,Rg_y,Rg_z,z[Nall-1]);

				// gyration_radius();
				
				// t_meas_temp = (float) t_meas*dt;
				// fwrite(&t_meas_temp, sizeof(t_meas_temp), 1, meas_bin_file);
				
				// count_sig_temp = (float) count_sig;
				// fwrite(&count_sig_temp, sizeof(count_sig_temp), 1, meas_bin_file);
				
				// count_en_temp = (float) count_en;
				// fwrite(&count_en_temp, sizeof(count_en_temp), 1, meas_bin_file);
				
				// Rg_x_temp = (float) Rg_x;
				// fwrite(&Rg_x_temp, sizeof(Rg_x_temp), 1, meas_bin_file);
				
				// Rg_y_temp = (float) Rg_y;
				// fwrite(&Rg_y_temp, sizeof(Rg_y_temp), 1, meas_bin_file);
				
				// Rg_z_temp = (float) Rg_z;
				// fwrite(&Rg_z_temp, sizeof(Rg_z_temp), 1, meas_bin_file);
				
				// z_end_temp = (float) z[Nall-1];
				// fwrite(&z_end_temp, sizeof(z_end_temp), 1, meas_bin_file);
				
				//fflush(meas_bin_file);
			}
		} /* t_meas */

		int lmnop;
		for(lmnop=0;lmnop<=(Nall-buffer)*Numb_meas;lmnop++){
			
			z_end_temp = (float)corrChain_z[lmnop];
			fwrite(&z_end_temp, sizeof(z_end_temp), 1, meas_bin_file);
			corrChain_z[lmnop] = 0.0;
			
			z_end_temp = (float)corrChain_nrg[lmnop];
			fwrite(&z_end_temp, sizeof(z_end_temp), 1, meas_bin_file_2);
			corrChain_nrg[lmnop] = 0.0;
		}

		finish();
		count = 0.0;

		t1 = clock(); dtime = (double) (t1 - t0);
		dtime /= CLOCK_PER_SEC;
		printf("Time:%f\n\n", dtime);
		//fclose(anim);
		
		fclose(meas_bin_file);
		fclose(meas_bin_file_2);

	}
	

	
	exit(0);
}
