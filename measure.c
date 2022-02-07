#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

double *rg, *g1, *g2, *g3, *g4, *g5, *cor;
double *x_end, *x_mid, *x_cm, *r_ee, *r2_ee, *l2_;
double *y_end, *y_mid, *y_cm, a2m;
double *z_end, *z_mid, *z_cm, am2, dnpm, dlpm;
double *e_kin, *e_pot, E_pot=0.0;
//Additions
double *rg0x, *rg0y, *rg0z;


long  put_measurment( char *fname )
{
long  i, it;
double rx, w, t, fac, Cv, w2;
FILE *fp;

dnpm = 1.0 / npm;
dlpm = 1.0 / lpm;

rx = dlpm*dnpm/nrun;
for(i=1, w=0.0; i<Numb_meas; i++) w +=rg[i];    rg[0] = w/(Numb_meas-1);
for(i=1, w=0.0; i<Numb_meas; i++) w +=l2_[i];   l2_[0] = w/(Numb_meas-1);
for(i=1, w=0.0; i<Numb_meas; i++) w +=r2_ee[i]; r2_ee[0] = w/(Numb_meas-1);

for(i=1, w=0.0, t=0.0, w2=0.0; i<Numb_meas; i++)
{
 w += e_kin[i];
 t += e_pot[i];
 w2 += square(e_kin[i]+e_pot[i])/(double)(Numb_meas-1);
}

e_kin[0] = w/(double)(Numb_meas-1);
e_pot[0] = t/(double)(Numb_meas-1);
Cv = (w2 - square(e_kin[0]+e_pot[0]))/(T*T)/Nall;

fp=fopen(fname,"w");
if(fp==(NULL)) return(-1);
//fprintf( fp,"acceptance: %lf Cv = %lf\n", accp/(nrun*mcsmax*Nmono),Cv);

am2 *= dnpm;
am2 /= (Numb_meas*nrun);
am2  = square(am2);
a2m *= (dnpm/(double)(Numb_meas*nrun));
   fac = dnpm /(double)(nrun*lpm1);

for( i=0; i<Numb_meas-1; i++){
   t  = 1.0/(double)((nrun*Numb_meas)-i-1);
   w  = dnpm*t;
   cor[i] = (cor[i]*t*dnpm);
   if((a2m-am2)!=0.0) cor[i] = ((cor[i]-am2)/(a2m-am2));
   t = (i+1)*step_meas;
   fprintf(fp,"%.3lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %f\n",
   t*dt, g1[i]*w, g3[i]*w, g4[i]*w, e_kin[i]*fac, e_pot[i]*fac, cor[i], rg[i]*rx, r2_ee[i]*fac*lpm1, l2_[i]*fac, r_ee[i]);
   }
fclose(fp);  
return(0);
}

void step_measure( long itime )
{
  long i, j, ic, k , it, jt, dt, i0, j0, mid1, mid2, m;
  double cmx, cmy, cmz, rg0, e, xee, yee, zee;
  double a12, dummy, ene, E_kin=0.0, half=0.5;


  dlpm=1.0/lpm; // number of beads in the chain


  it = itime % Numb_meas;

  mid1 = lpm/2;
  mid2 = mid1 + 1;
 
for(ic=0; ic<npm; ic++){
   i = ic * lpm;
   k = (it * npm) + ic; // npm = # polymer chains

   x_end[k]=x[i]; y_end[k] = y[i]; z_end[k] = z[i];

   x_mid[k]=(x[i+mid1]+x[i+mid2])/2;
   y_mid[k]=(y[i+mid1]+y[i+mid2])/2;
   z_mid[k]=(z[i+mid1]+z[i+mid2])/2;

   xee = x[i] - x[i+lpm1];
   yee = y[i] - y[i+lpm1];
   zee = z[i] - z[i+lpm1];

   rg0 = square(xee) + square(yee) + square(zee);
   r_ee[k] =  sqrt(rg0);
   if(itime % Numb_meas)   r2_ee[it] += rg0;

   a2m += square( r_ee[k] );
   am2 += r_ee[k];
   cmx=cmy=cmz=0.0;

   for(j=0; j<lpm; j++){
      cmx += x[i+j]; cmy += y[i+j]; cmz += z[i+j];
      }
   cmx *= dlpm; cmy *= dlpm; cmz *= dlpm;
   x_cm[k] = cmx;
   y_cm[k] = cmy;
   z_cm[k] = cmz;

   if(itime % Numb_meas) {
      rg0 = 0.0;
      for(j=0; j<lpm; j++){
         rg0 += square(x[i+j]-cmx);
         rg0 += square(y[i+j]-cmy);
         rg0 += square(z[i+j]-cmz);
      }
      rg[it] += rg0;

  // calculate kin energy
  E_kin = 0.0;
  for (m=0;m<lpm1;m++) {
    ene = square(vx[i+m]) + square(vy[i+m]) + square(vz[i+m]);
    E_kin += ene*half*mass;
  }
  e_kin[it] += E_kin;

  // compute kinetic energy
  calculate_potential_energy(ic);
  e_pot[it] += E_pot;

      for(j=0; j<lpm1; j++)
         l2_[it] += square(x[i+j]-x[i+j+1])+square(y[i+j]-y[i+j+1])+square(z[i+j]-z[i+j+1]);
     }
   } // ic   
    
for(dt=1; dt<min(Numb_meas,itime); dt++){
   jt = dt - 1;
   a12 = 0.0;
   j0 = (it*npm);
   i0 = ((Numb_meas + it - dt) % Numb_meas)*npm;    //(((it-dt))*npm); 
   for(ic=0; ic<npm; ic++){
      j = j0 + ic; i = i0 + ic;
      g1[jt] += square(x_mid[j]-x_mid[i]);
      g1[jt] += square(y_mid[j]-y_mid[i]);
      g1[jt] += square(z_mid[j]-z_mid[i]);
//      g2[jt] += square(x_mid[j]-x_cm[j]+x_cm[i]-x_mid[i]);
//      g2[jt] += square(y_mid[j]-y_cm[j]+y_cm[i]-y_mid[i]);
//      g2[jt] += square(z_mid[j]-z_cm[j]+z_cm[i]-z_mid[i]);
      g3[jt] += square(x_cm[j]-x_cm[i]);
      g3[jt] += square(y_cm[j]-y_cm[i]);
      g3[jt] += square(z_cm[j]-z_cm[i]);
      g4[jt] += square(x_end[j]-x_end[i]);
      g4[jt] += square(y_end[j]-y_end[i]);
      g4[jt] += square(z_end[j]-z_end[i]); 
//      g5[jt] += square(x_end[j]-x_cm[j]+x_cm[i]-x_end[i]);
//      g5[jt] += square(y_end[j]-y_cm[j]+y_cm[i]-y_end[i]);
//      g5[jt] += square(z_end[j]-z_cm[j]+z_cm[i]-z_end[i]);

      a12 += r_ee[j]*r_ee[i];
   } 
   cor[jt] += a12;
 }
} 


void init_measure(void)
{
long itime;

itime = Numb_meas+1;

rg = (double *) calloc( itime, sizeof(double) );
r2_ee = (double *) calloc( itime, sizeof(double) );
l2_ = (double *) calloc( itime, sizeof(double) );

e_kin = (double *) calloc( (size_t) itime, sizeof(double) );
e_pot = (double *) calloc( (size_t) itime, sizeof(double) );

g1 = (double *) calloc( itime, sizeof(double) );
g2 = (double *) calloc( itime, sizeof(double) );
g3 = (double *) calloc( itime, sizeof(double) );
g4 = (double *) calloc( itime, sizeof(double) );
g5 = (double *) calloc( itime, sizeof(double) );

cor    = (double *) calloc( itime, sizeof(double) );

itime = Numb_meas*npm;

x_end = (double *) calloc( itime, sizeof(double) );
y_end = (double *) calloc( itime, sizeof(double) );
z_end = (double *) calloc( itime, sizeof(double) );
x_mid = (double *) calloc( itime, sizeof(double) );
y_mid = (double *) calloc( itime, sizeof(double) );
z_mid = (double *) calloc( itime, sizeof(double) );
r_ee  = (double *) calloc( itime, sizeof(double) );
x_cm = (double *) calloc( itime, sizeof(double) );
y_cm = (double *) calloc( itime, sizeof(double) );
z_cm = (double *) calloc( itime, sizeof(double) );

if(rg==NULL) goto error;
if(r2_ee==NULL) goto error;
if(l2_==NULL) goto error;

if(g1==NULL) goto error;
if(g2==NULL) goto error;
if(g3==NULL) goto error;
if(g4==NULL) goto error;
if(g5==NULL) goto error;

if(cor==NULL) goto error;

if(x_end==NULL) goto error;
if(y_end==NULL) goto error;
if(z_end==NULL) goto error;
if(x_mid==NULL) goto error;
if(y_mid==NULL) goto error;
if(z_mid==NULL) goto error;
if(r_ee==NULL) goto error;
if(x_cm==NULL) goto error;
if(y_cm==NULL) goto error;
if(z_cm==NULL) goto error;

goto quit;

error:
printf("*** ERROR during allocation of the Stat arrays\n");
exit(-1);

quit:
return;
}


void calculate_potential_energy(long ic)
{
  long i, j, k, ip;
  double dr2, ene, half = 0.5;

// Compute potential energy from atom positions
  E_pot = 0.0;
  // compute bonded & non-bonded interactions from the Verlet-list

  for(i=ic;i<ic+lpm;i++) {
    ip = i+1;
    for(k=Marker1[i];k<=Marker2[i];k++) { //search in Verlet-list
      j = List[k];
      dr2 = square(x[i]-x[j]) + square(y[i]-y[j]) + square(z[i]-z[j]);
//    ***************************************************
      if ((j==ip)&&(j%lpm)){ // right neighbor along the backbone!
        if (Morse) { // if NN are bonded by Morse potential
          ene = Morse_pot(dr2);
          E_pot += ene ;
          }
        else { // if NN are bonded by Kremer-Grest potential
          ene = KG_pot(dr2);
          E_pot += ene ;
        }
      } /* right NN */
      else {  // is j a non-bonded neighbor?
        if ((WCA)&&(dr2 <= Range_WCA2)) {
          ene = WCA_pot(dr2);
          E_pot += ene ;
        }
        if ((LJ)&&(dr2 <= Range_LJ2)) {
          ene = LJ_pot(dr2);
          E_pot += ene ;
        }
      } /* else j is non-bonded meighbor */
    } /* for k */

  } /* for i */
  return;
}
