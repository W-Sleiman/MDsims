#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

#define square(x) ((x)*(x))

long get_params(char *fname)
{ // read main parameters from the input file "in"
 FILE *fp;
 char buf[120];

 fp=fopen(fname,"r");
 if(fp==(NULL)) return(-1);
 long seed1;
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&nrun);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&npm);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&lpm);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&steps);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&step_meas);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&dt);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&mass);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&T);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&Gamma);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&seed1);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&Morse);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&e_M);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&alpha);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&r_min);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&LJ);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&WCA);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&e_LJ);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&sigma);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&KG);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&k_spring);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&R0);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&f_ext);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&skin);
 fgets(buf,120,fp);   sscanf(buf,"%s ", fin1);
 fgets(buf,120,fp);   sscanf(buf,"%s ", fout1);
 fgets(buf,120,fp);   sscanf(buf,"%s ", fout2);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&gen_init);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&e_b);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&sigma_b);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&time_width);
 fgets(buf,120,fp);   sscanf(buf,"%lf ",&nrg_ncrmnt);
 fgets(buf,120,fp);   sscanf(buf,"%ld ",&buffer);
 
 fclose(fp);
 if(nrg_step == 0){
	seed = -need;
 }
 //seed = (long) -3;
 return(0);
}

long  get_monomers(char *fname)
{ // read monomer positions and velocities 
 long  i;
 double rx, ry, rz, Vx, Vy, Vz;
 FILE *fp;
 char buf[80];
 if(!gen_init){
	fp=fopen(fname,"r");
	if(fp==(NULL)) return(-1);
	for(i=0; i< Nall; i++){
		fgets(buf,80,fp);
		sscanf(buf,"%lf,%lf,%lf,%lf,%lf,%lf", &rx, &ry, &rz, &Vx, &Vy, &Vz);
		x[i] = rx; y[i]=ry; z[i]=rz; vx[i] = Vx; vy[i]=Vy; vz[i]=Vz;
	}
	fclose(fp);
	return(0);
 } 
 else if(gen_init){
	 
	fp=fopen(fname,"w");
	if(fp==(NULL)) return(-1);

	for(i=0; i< Nall; i++){
		x[i]=0.0; y[i]=0.0; z[i]=i*sigma; 
		vx[i] = ran1(&seed)-0.5; vy[i]= ran1(&seed)-0.5; vz[i]= ran1(&seed)-0.5;
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",x[i],y[i],z[i],vx[i],vy[i],vz[i]);
		
	}
	/*for(i=4; i< Nall; i++){
		x[i]=0.0; y[i]=(i-3)*sigma; z[i]=z[i-1]; 
		vx[i] = ran1(&seed)-0.5; vy[i]= ran1(&seed)-0.5; vz[i]= ran1(&seed)-0.5;
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",x[i],y[i],z[i],vx[i],vy[i],vz[i]);
		
	}*/
	fclose(fp);
	return(0);
 }
}

long  put_monomers(char *fname)
{ // store monomer final positions and final velocities
 long  i;
 double Rx, Ry, Rz, Vx, Vy, Vz;
 FILE *fp;

 fp=fopen(fname,"w");
 if(fp==(NULL)) return(-1);

 
 for(i=0; i<Nall; i++){
   Rx=x[i]; Ry=y[i]; Rz=z[i]; Vx=vx[i]; Vy=vy[i]; Vz=vz[i];
   fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n", Rx, Ry, Rz, Vx, Vy, Vz);
 }
 fclose(fp);
 return(0);
}

long init_mem()

  // allocate memory for arrays
{

 f_aux = (double *) calloc((size_t)  Nall, sizeof(double) );
 f_cons = (double *) calloc((size_t)  Nall, sizeof(double) );

 x = (double *) calloc((size_t) Nall, sizeof(double) );
 y = (double *) calloc((size_t) Nall, sizeof(double) );
 z = (double *) calloc((size_t) Nall, sizeof(double) );

 xold = (double *) calloc((size_t) Nall, sizeof(double) );
 yold = (double *) calloc((size_t) Nall, sizeof(double) );
 zold = (double *) calloc((size_t) Nall, sizeof(double) );

 vx = (double *) calloc( (size_t) Nall, sizeof(double) );
 vy = (double *) calloc((size_t)  Nall, sizeof(double) );
 vz = (double *) calloc((size_t)  Nall, sizeof(double) );


 cons_fx = (double *) calloc( (size_t) Nall, sizeof(double) );
 cons_fy = (double *) calloc((size_t)  Nall, sizeof(double) );
 cons_fz = (double *) calloc((size_t)  Nall, sizeof(double) );

 fx_aux = (double *) calloc( (size_t) Nall, sizeof(double) );
 fy_aux = (double *) calloc((size_t)  Nall, sizeof(double) );
 fz_aux = (double *) calloc((size_t)  Nall, sizeof(double) );

 List = (long *) calloc((size_t) MaxPairsPerAtom*Nall, sizeof(long) );

 Marker1 = (long *) calloc((size_t)  Nall, sizeof(long) );
 Marker2 = (long *) calloc((size_t)  Nall, sizeof(long) );
 Advance = (long *) calloc((size_t)  Nall, sizeof(long) );
 
 DisplaceList_x = (double *) calloc((size_t)  Nall, sizeof(double) );
 DisplaceList_y = (double *) calloc((size_t)  Nall, sizeof(double) );
 DisplaceList_z = (double *) calloc((size_t)  Nall, sizeof(double) );
 //NEW
 index_ads  = (double *) calloc((size_t) Nall-1, sizeof(double) );
 index_ads_energy  = (double *) calloc((size_t) Nall-1, sizeof(double) );
 Zcor  = (double *) calloc((size_t)Numb_meas+1, sizeof(double));
 Zcor2 = (double *) calloc((size_t)Numb_meas+1, sizeof(double));
 Mcor  = (double *) calloc((size_t)Numb_meas+1, sizeof(double));
 Mcor2 = (double *) calloc((size_t)Numb_meas+1, sizeof(double));
 
 corrChain_z = (double *) calloc((size_t)(Nall-buffer)*Numb_meas, sizeof(double));
 corrChain_nrg = (double *) calloc((size_t)(Nall-buffer)*Numb_meas, sizeof(double));

 // ==============================================================

 if( x   ==NULL ) goto error;
 if( y   ==NULL ) goto error;
 if( z   ==NULL ) goto error;

 if( xold   ==NULL ) goto error;
 if( yold   ==NULL ) goto error;
 if( zold   ==NULL ) goto error;

 if( vx   ==NULL ) goto error;
 if( vy   ==NULL ) goto error;
 if( vz   ==NULL ) goto error;

 if( cons_fx==NULL ) goto error;
 if( cons_fy==NULL ) goto error;
 if( cons_fz==NULL ) goto error;

 if( fx_aux==NULL ) goto error;
 if( fy_aux==NULL ) goto error;
 if( fz_aux==NULL ) goto error;

 if( List   ==NULL ) goto error;
 if( Marker1==NULL ) goto error; 
 if( Marker2==NULL ) goto error; 
 if( Advance==NULL ) goto error; 

 if( DisplaceList_x==NULL ) goto error;
 if( DisplaceList_y==NULL ) goto error;
 if( DisplaceList_z==NULL ) goto error;

//NEW
 if( index_ads==NULL ) goto error;
 if( index_ads_energy==NULL ) goto error;
 if( Zcor ==NULL ) goto error;
 if( Zcor2==NULL ) goto error;
 if( Mcor ==NULL ) goto error;
 if( Mcor2==NULL ) goto error;
 if( corrChain_z==NULL ) goto error;
 if( corrChain_nrg==NULL ) goto error;
//
 return(0);
 error:
   puts("*** ERROR while allocating arrays");
   exit(-1);
   return(0);
}

void initialize(void)
{
 long i, errflg;
 double ksi;
 // read main parameters from input file "in" +++++++++++++++++++++ 
 errflg  =  get_params( fcfg ); 
 if(errflg) { puts("*** ERROR in 'in' file INPUT operations"); exit(-1);}
 
 //printf("1: e_b = %lg \n",e_b);
 
 if (Morse && KG) { puts("*** ERROR: A bond must be either Morse OR KG"); 
   exit(-1);}
 else if (WCA && LJ) { puts("*** ERROR: A non-bond must be either WCA OR LJ");
   exit(-1);}
 
 //printf("2: npm = %ld,	lpm = %ld \n",npm, lpm);
 
 e_b = e_b + nrg_ncrmnt * nrg_step; 
 Nall = npm * lpm; // total number of beads
 lpm1 = lpm -1;
 Numb_meas = steps*nrun / step_meas;

 //printf("3: npm = %ld,	lpm = %ld \n",npm, lpm);
 //printf("4: Nall = %ld \n",Nall);
 if (e_b>2){puts("*** ERROR: maximum e_b value (which is 2.0) exceeded ***"); exit(-1);}


 // allocate memory for arrays ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 init_mem();
 
 //printf("5: Nall = %ld \n",Nall);
 
 ListUpdateRequested = 1; // update Verlet_list in the beginning
 MaxListLength = MaxPairsPerAtom * Nall; // fix the Verlet_list length
 
 // if(e_b<0.15){	
 // sprintf(fin1,"initial_conditions/P_end_%lg_%ld",0.01,50);
 // sprintf(fin1,"initial_conditions/P_end_0_50");
 // }
 // else if (e_b>=0.15){
 sprintf(fin1,"initial_conditions/P_end_%lg_%ld",e_b,Nall);
 // }
 
 // read starting configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 errflg |= get_monomers(fin1);
 if(errflg) { puts("*** ERROR in 'IC' file INPUT operations"); exit(-1);}
 
 //printf("6: z(1) = %lg \n",z[1]);
 //printf("7: vz(1) = %lg \n",vz[1]);


 // compute necessary constants for the Langevin thermostat

 inv_mass = 1.0 / mass;
 inv_sigma2 = 1.0 / (sigma * sigma);
 sqr_Gamma_dt = sqrt(24 * T * Gamma * mass / dt);

 Range_WCA = pow(2.0,(1.0/6.0));
 Range_WCA2 = Range_WCA * Range_WCA;
 Range_LJ = 2.5 * sigma;
 Range_LJ2 = Range_LJ * Range_LJ;
 
 inv_sigma_b2 = 1.0 / (sigma_b * sigma_b);
 Range_LJ_b = 2.5 * sigma_b;
 
 LJ_shift = 0.25;
 if (WCA) {
   Range_Verlet = Range_WCA + skin;
 }
 else {
   Range_Verlet = 2.5*sigma + skin;
 }
}

void finish(void)
{
 long errflg;

 // store configuration + results

 errflg = put_monomers(fout1);
 //errflg |= put_measurment(fout2);
 if(errflg) { puts("*** ERROR in file OUTPUT operations"); exit(-1);}
 puts("*** That's all folks...");
}

void f_list(void)
{ // compute forces acting on monomer i
  long i;
  for(i=0;i<Nall;i++) {
    f_cons[i] = square(cons_fx[i])+ square(cons_fy[i])+ square(cons_fz[i]);
    f_aux[i] = square(fx_aux[i])+ square(fy_aux[i])+ square(fz_aux[i]);
  }
}

void check_bonds()
{  // compute squared bond lengths
  long i, k;
  double bond2;
  
  for(i=0;i<Nall;i++) {
    if (i%lpm1) {
      bond2 = square(x[i+1]-x[i])+square(y[i+1]-y[i])+square(z[i+1]-z[i]);
      if(bond2 >= R0*R0) {
        printf("error: bond %ld = %lf too long! at time: %ld\n",i,sqrt(bond2),t_meas);
		finish();
		exit(-1);
      }
	  if(isnan(x[i])){
		printf("error: nan x %ld at time: %ld\n",i,t_meas);
		finish();
		exit(-1);
	  }
    } // i must not be the last bead in a chain
  } /* i */
}

void cons_force(void)
{
  // Calculate conservative forces
  long i, j, k, ip;
  double dr2, isotrop_f, fcx, fcy, fcz;

  for(i=1;i<Nall;i++) {
    cons_fx[i] = cons_fy[i] = cons_fz[i] = 0.0;
	double dz2 = z[i]*z[i];
	if (z[i] <= Range_LJ_b){ // 
        cons_fz[i] += force_LJ_boundary(dz2)*z[i];
      } 
  }

  for(i=0;i<Nall-1;i++) {
    ip = i+1;
    for(k=Marker1[i];k<=Marker2[i];k++) { //search in Verlet-list
      j = List[k];
      dr2 = square(x[i]-x[j]) + square(y[i]-y[j]) + square(z[i]-z[j]);
	  
      fcx = fcy = fcz = 0.0;
      
//   *******************************************************
      if ((j==ip)&&(j%lpm)){ // if right NN along the backbone!
        if (Morse) { // if NN is bonded by Morse potential
          isotrop_f = force_M(dr2);
            fcx += isotrop_f * (x[i] - x[j]);
            fcy += isotrop_f * (y[i] - y[j]);
            fcz += isotrop_f * (z[i] - z[j]);
        } /* if Morse */
        else if(KG) { // if NN is bonded by Kremer-Grest potential
          isotrop_f = force_KG(dr2);
            fcx += isotrop_f * (x[i] - x[j]);
            fcy += isotrop_f * (y[i] - y[j]);
            fcz += isotrop_f * (z[i] - z[j]);
        } /* if KG */
      } /* if right NN */

      else {  // j is not a NN
        if (WCA) {
          if (dr2 <= Range_WCA2) {
            isotrop_f = force_WCA(dr2);
            fcx += isotrop_f * (x[i] - x[j]);
            fcy += isotrop_f * (y[i] - y[j]);
            fcz += isotrop_f * (z[i] - z[j]);
          }
        } /* if WCA */  
        if (LJ) {
          if (dr2 <= Range_LJ2) {
            isotrop_f = force_LJ(dr2);
            fcx += isotrop_f * (x[i] - x[j]);
            fcy += isotrop_f * (y[i] - y[j]);
            fcz += isotrop_f * (z[i] - z[j]);
          }
        } /* if LJ */ 
      } /* j is not NN */
          
      cons_fx[i]+=fcx; cons_fy[i]+=fcy; cons_fz[i]+=fcz;
      cons_fx[j]-=fcx; cons_fy[j]-=fcy; cons_fz[j]-=fcz;

    } /* for k */
  } /* for i */
}

//======================================================
double force_M(double x2)
 // Compute Morse force
{
 double two=2.0, dr, f, F_M;

 dr = sqrt(x2);
 f = exp(-alpha * (dr - r_min));
 F_M = two * alpha * e_M * (f*f - f) / dr;
 
 return(F_M);
}
//======================================================
double force_KG(double x2)
 // Compute Kremer-Grest force: x2 is squared distance!
{
 double two=2.0, f, f_wca=0.0, s_to_x2, s_to_x6, F_KG;

 f = - k_spring * R0*R0/(R0*R0 - x2);
 F_KG = f;
 if (x2 <= Range_WCA2) {
   s_to_x2 = sigma * sigma / x2;
   s_to_x6 = s_to_x2 * s_to_x2 * s_to_x2;
   f_wca = 24 * e_LJ * (two * s_to_x6 * s_to_x6 - s_to_x6)*s_to_x2; 
   f_wca *= inv_sigma2;
 }
 F_KG += f_wca;
 return(F_KG);
}

//======================================================
double force_WCA(double x2)
 // Compute Lennard-Jones force: x2 is squared distance!
 {
  double two=2.0, f, s_to_x2, s_to_x6, F_WCA=0.0;
  s_to_x2 = sigma * sigma / x2;
  s_to_x6 = s_to_x2 * s_to_x2 * s_to_x2;
  f = 24 * e_LJ * (two * s_to_x6 * s_to_x6 - s_to_x6)*s_to_x2;
  F_WCA = f*inv_sigma2;
  return(F_WCA);
 }
       
//======================================================
double force_LJ(double x2)
 // Compute Lennard-Jones force: x2 is squared distance!
{
 double two=2.0, f, s_to_x2, s_to_x6, F_LJ;

 s_to_x2 = sigma * sigma / x2;
 s_to_x6 = s_to_x2 * s_to_x2 * s_to_x2;
 f = 24 * e_LJ * (two * s_to_x6 * s_to_x6 - s_to_x6)*s_to_x2; 
 F_LJ = f*inv_sigma2;
 return(F_LJ);
}
//======================================================
double force_LJ_boundary(double x2)
 // Compute Lennard-Jones force: x2 is squared distance!
{
 double two=2.0, f_b, s_to_x2, s_to_x6, F_LJ_b;

 s_to_x2 = sigma_b * sigma_b / x2;
 s_to_x6 = s_to_x2 * s_to_x2 * s_to_x2;
 f_b = 24 * e_b * (two * s_to_x6 * s_to_x6 - s_to_x6)*s_to_x2; 
 F_LJ_b = f_b*inv_sigma_b2;
 return(F_LJ_b);
}
//======================================================
double Morse_pot(double x2)
 // here x2 is squared distance!
 // Compute Morse potential energy for r_min, e_min and alpha known! 
 // At dr = r_min U_M = 0 !!! (and not - e_M)
{
  double one = 1.0, dr, ex, U_M;
  dr = sqrt(x2);
  
  ex = one - exp(-alpha * (dr - r_min));
  U_M = e_M * ex * ex;
  return(U_M); 
}
//======================================================
double KG_pot(double x2)
 // Compute Kremer-Grest potential (FENE + WCA)
 // FENE elastic constant is k_spring, range is R0
{ 
  double one = 1.0, half = 0.5, four = 4.0, U_FENE, U_WCA=0.0;
  double s_to_x, LJ_shift = 0.25, u;
 
  U_FENE = - half * k_spring * R0 * R0 * log(one - x2/(R0*R0));
  if (x2 <= Range_WCA2) {  
    s_to_x = sigma*sigma/x2; // here x2 is the distance squared!
    s_to_x = s_to_x * s_to_x * s_to_x;
    U_WCA = four * e_LJ * ((s_to_x * s_to_x - s_to_x) + LJ_shift);
  }
  u = U_FENE + U_WCA;
  return(u);
}
//======================================================
double LJ_pot(double x2)
 // Compute Lennard-Jones potential: x2 is distance squared! 
{
 double four = 4.0, s_to_x, U_LJ; 

 s_to_x = sigma*sigma/x2;
 s_to_x = s_to_x * s_to_x * s_to_x;
 U_LJ = four * e_LJ * (s_to_x * s_to_x - s_to_x);

 return(U_LJ);
}
//======================================================
double WCA_pot(double x2)
 // Compute Weeks-Chandler-Anderson potential: x = distance^2
{
 double four = 4.0, s_to_x, U_WCA=0.0;
 if (x2 <= Range_WCA2) {
   s_to_x = sigma*sigma/x2;
   s_to_x = s_to_x * s_to_x * s_to_x;
   U_WCA = four * e_LJ * ((s_to_x * s_to_x - s_to_x) + LJ_shift);
 }
 return(U_WCA);
}
//======================================================
double LJ_pot_boundary(double x2)
 // Compute Lennard-Jones potential: x2 is distance squared! 
{
 double four = 4.0, s_to_x, U_LJ_b; 

 s_to_x = sigma_b*sigma_b/x2;
 s_to_x = s_to_x * s_to_x * s_to_x;
 U_LJ_b = four * e_b * (s_to_x * s_to_x - s_to_x);

 return(U_LJ_b);
}

void get_gasdev(long iseed)
{
  long i;
  double rrr=1.0;

  for(i=0; i< (2*Nall); i++){
    gauss[i]=0.0;
    gauss[i]=alpha*gasdev(&seed);
  }
}

/* void adsorption_measurement(){
	//Still need to declare time_width and tell vvmd how to call this method
	int i;
	double ad_cnt = 0.0, EE = 0.0;
	for(i=1;i<Nall;i++){
		if(z[i]<Range_LJ_b){
			ad_cnt = ad_cnt + 1.0;
			//EE = EE + LJ_pot_boundary(z[i]*z[i]);			
		}
	}
	//count2 += ad_cnt*ad_cnt;
	count = ad_cnt;
	//E += EE/(2.0*time_width*ad_cnt);
	
} */

void adsorption_measurement(){
	//Still need to declare time_width and tell vvmd how to call this method
	int i;
	double ad_cnt = 0.0, EE = 0.0;
	for(i=1;i<Nall;i++){
		if(z[i]<Range_LJ_b){
			//ad_cnt = ad_cnt + 1.0;
			EE = EE + LJ_pot_boundary(z[i]*z[i]);			
		}
	}
	//count2 += ad_cnt*ad_cnt;
	count = EE/e_b;
	//E += EE/(2.0*time_width*ad_cnt);
	
}

void gyration_radius(){
	int ii, jj;
	double zcm = 0.0, ycm = 0.0, xcm = 0.0; 
	double rgx = 0.0, rgy = 0.0, rgz = 0.0, ad_cnt = 0.0, EE = 0.0; 
	
	for(ii=1;ii<Nall;ii++){
		//zcm += z[ii]/Nall;
		ycm += y[ii]/Nall;
		xcm += x[ii]/Nall;
		zcm += z[ii]/Nall;
		if(z[ii]<Range_LJ_b){
			ad_cnt = ad_cnt + 1.0;
			EE = EE + LJ_pot_boundary(z[ii]*z[ii]);			
		}
	}
	count_sig = ad_cnt;
	count_en = EE/e_b;
	for(jj=1;jj<Nall;jj++){
		rgx += fabs(x[jj]-(xcm))/Nall;
		rgy += fabs(y[jj]-(ycm))/Nall;
		rgz += fabs(z[jj]-(zcm))/Nall;
		//printf("hello monomer %d : rgz is %lg, z[jj] = %lg, zcm = %lg   at time %d ;\n",jj,rgz,z[jj],zcm,t_meas);
		//rgxy += sqrt(square(x[jj]-xcm)+square(y[jj]-ycm))/Nall;
	}
	//Rg_perp += rgz /(2.0*time_width);
	Rg_x  = rgx;
	Rg_y  = rgy;
	Rg_z  = rgz;

}

void Corr_Along_Chain(){
	//Need to define and alloc corrChain (array of correlations along chain) and buffer (#of monomers skipped near grafted end)
	int i;
	int j;
	
	for(j=0;j<Nall-buffer;j++){
		//for(j=0;j<Nall-buffer-i;j++){
		if(z[j+buffer]<Range_LJ_b){
			corrChain_z [j+((t_meas/step_meas)-1)*(Nall-buffer)] = corrChain_z[j+((t_meas/step_meas)-1)*(Nall-buffer)] + 1;
		}
		corrChain_nrg [j+((t_meas/step_meas)-1)*(Nall-buffer)] = corrChain_nrg[j+((t_meas/step_meas)-1)*(Nall-buffer)] + LJ_pot_boundary(z[j+buffer]*z[j+buffer])/e_b;
		
			// if(z[i+buffer]<Range_LJ_b && z[i+j+buffer]<Range_LJ_b){
				// corrChain[j] = corrChain[j] + 1;
				// printf("1: z[i+buffer] = %lg,	z[i+j+buffer] = %lg \n",z[i+buffer], z[i+j+buffer]);
				// printf("2: corrChain[j] = %lg,	\n",corrChain[j]);
				// /(Numb_meas*nrun)
				// /(e_b*e_b)
			//}
		
		//}
	}
	
}
