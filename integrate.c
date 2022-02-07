#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

void velocity_Verlet()
{ //Langevin Dynamics velocity-Verlet algorithm

  long j;
  double random, noise, half = 0.5;

  if(ListUpdateRequested) {
    Verlet_List();
    ListUpdateRequested = 0;
  }

  // Do 1-st half integration step
  //==========================
 
  for(j=1;j<Nall;j++) {

    xold[j] = x[j]; yold[j] = y[j]; zold[j] = z[j];

  // 1. update velocities at 1/2-dt
  //v(t+dt/2) = v(t) + f(t)*dt/(2m)
    vx[j] += half * inv_mass * fx_aux[j] * dt;
    vy[j] += half * inv_mass * fy_aux[j] * dt;
    vz[j] += half * inv_mass * fz_aux[j] * dt;

  // 2. full update of monomer coordinates
  // x(t+dt) = x(t) + v(t+dt/2)*dt
    x[j] += vx[j]*dt;
    y[j] += vy[j]*dt;
    z[j] += vz[j]*dt;
  } 
  
  /*if(!(t_meas%1000)){
	  for (j = 0; j<Nall; j++){
		  fprintf(anim, "%ld,%lf,%lf,%lf\n", t_meas, x[j], y[j], z[j]);
	  }
  }*/
  
  // 3. Compute forces at the new x, y, z
  // f(t+dt) = f_cons(t) - Gamma*v(t+dt/2) + sqrt(24*T*Gamma*m/dt)*RAND
    cons_force();
  
  // 4. update velocities at full step dt

    for(j=1;j<Nall;j++) {

      random = ran1(&seed);
      noise = sqr_Gamma_dt * (random-half);
      fx_aux[j] = cons_fx[j] - mass*Gamma*vx[j] + noise;

      random = ran1(&seed);
      noise = sqr_Gamma_dt * (random-half);
      fy_aux[j] = cons_fy[j] - mass*Gamma*vy[j] + noise;

      random = ran1(&seed);
      noise = sqr_Gamma_dt * (random-half);
      fz_aux[j] = cons_fz[j] - mass*Gamma*vz[j] + noise;

      vx[j] += half * inv_mass * fx_aux[j] * dt;
      vy[j] += half * inv_mass * fy_aux[j] * dt;
      vz[j] += half * inv_mass * fz_aux[j] * dt;
    } /* end of integration over dt */ 
     
  //=======================================
  // Update array DisplaceList which contains 
  // atoms displacements since last update

  for(j=1;j<Nall;j++) {
    DisplaceList_x[j] += x[j]-xold[j];
    DisplaceList_y[j] += y[j]-yold[j];
    DisplaceList_z[j] += z[j]-zold[j]; 
  }
  // if nececcery, update Verlet_List
  ListUpdateRequested = MovedTooMuch();
}

void Verlet_List(void)
{
  long i, j, L=0;
  double dr2, Range_Verlet_sq;
  // Fincham-Ralston loop for list updating
  Range_Verlet_sq = Range_Verlet * Range_Verlet;
  
  for (i=0;i<Nall-1;i++) {
    for(j=i+1;j<Nall;j++) {

      dr2 = square(x[i]-x[j]) + square(y[i]-y[j]) + square(z[i]-z[j]);

      // if j is a neighbor, put 1, otherwise 0
      if (dr2 <= Range_Verlet_sq) {Advance[j] = 1;}
      else {Advance[j] = 0;}
    }
    Marker1[i] = L; // List for i starts here

    for(j=i+1;j<Nall;j++) {
      if (L > MaxListLength) goto error;
      List[L] = j;
      L += Advance[j];
    }
    Marker2[i] = L-1;

  } /* end i */
  for (i=0;i<Nall;i++) {
    // initialize monomer displacements
    DisplaceList_x[i] = 0.0; 
    DisplaceList_y[i] = 0.0; 
    DisplaceList_z[i] = 0.0; 
  }
  return;

  error:
  puts("*** ERROR: Verlet_list too small!");
  exit(-1);
  return;
}

int MovedTooMuch()
{ // compute maximum displacement of an atom 
  long i=0, MTM=0;
  double disp_1=0.0,disp_2=0.0, disp;
  
  for(i=1;i<Nall;i++) {
    disp = square(DisplaceList_x[i])+square(DisplaceList_y[i])+square(DisplaceList_z[i]);
    disp = sqrt(disp);
    if (disp >= disp_1) {
      disp_2 = disp_1;
      disp_1 = disp;
    }
    else if (disp_1 >= disp_2) disp_2 = disp_1;
  }         
  if (disp_1+disp_2 > skin) {
    MTM=1; // Verlet-list needs update
  }

  return(MTM);
}
