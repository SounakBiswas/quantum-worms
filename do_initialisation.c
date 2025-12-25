#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#define  next_nbr(ii,jj,kk,xx,yy,zz)  ((xx+ii+lx)%lx + ((yy+jj+ly)%ly)*lx + ((zz+lz+kk)%lz)*lx*ly)
void do_initialisation(void)
{
  double genrand_real2(void);
  //void init_by_array( unsigned long *init_key, unsigned long key_length);
  void init_genrand(unsigned long s);
  double choice;
  unsigned int i,j,k,site,nextsite;
  //unsigned long init[4]=  {0x772, 0x456, 0x110, 0x864}  /* {0x446,0x322,0x768,0x420}*/ , length2=4;
  unsigned long seed;
  opstr_l=INIT_OPSTR_L;
  opstr=(int*)malloc(opstr_l*sizeof(int));
  divider=(int*)malloc(opstr_l*sizeof(int));
  minority=(int*)malloc(opstr_l*sizeof(int));
  //init_by_array(init, length2);
  for(i=0;i<opstr_l;i++)
  opstr[i]=divider[i]=minority[i]=-1;
  n_niop=0;
  n_triagop=0;
  n_flipop=0;
  seed=SEED;
  init_genrand(seed);
  take_stat_bool=0;
  num_cluster=0;
  if(initialstate==0){ //infinite temperature initialization
    for(i=0;i<lx;i++){
      for(j=0;j<ly;j++){
        site=i+lx*j;
        choice=genrand_real2();
        if(choice < 0.5){
          sigma[site]=-1;
        } else{ sigma[site]=1;}
      }
    }


  }//initialstate=0

  if(initialstate==1){ //zero-temp initialization. (lx,ly mult. of 3)
    for(i=0;i<lx;i++){
      for(j=0;j<ly;j++){
        site=i+lx*j;
        if(((i+j)%3) == 0 ){
          sigma[site]=+1;  
        } else{
          sigma[site]=-1;
        }
      }
    }


  }//initialstate=1

  if(initialstate==2){ //zero temp (lx,ly even)
    for(i=0;i<lx;i++){
      for(j=0;j<ly;j++){
        site=i+lx*j;
        if((j%2) ==0 ){
          sigma[site]=+1;  
        } else{
          sigma[site]=-1;
        }
      }
    }


  }//initialstate=2


  if(initialstate==4){ //zero temp (lx,ly even)
    for(i=0;i<lx;i++){
      for(j=0;j<ly;j++){
        site=i+lx*j;
        if((i%2) ==0 ){
          sigma[site]=+1;  
        } else{
          sigma[site]=-1;
        }
      }
    }


  }//initialstate=4

  if(initialstate==6){ //zero temp (lx,ly even)
    for(i=0;i<lx;i++){
      for(j=0;j<ly;j++){
        site=i+lx*j;
        if(((i-j+ly)%2) ==0 ){
          sigma[site]=+1;  
        } else{
          sigma[site]=-1;
        }
      }
    }


  }//initialstate=4

//global measurements
  energy_av=0.0;
  sigmax_q_stat_av=sigmaz_q_stat_av=sigmaz_q_et_av=sigmaz_0_et_av=sigmaz_0_stat_av=0.0;
//local measurements

tautau_et_lby2_av=0.0;tautau_et_lby3_av=0.0;tautau_et_lby6_av=0.0;tautau_stat_lby2_av=0.0;tautau_stat_lby3_av=0.0;tautau_stat_lby6_av=0.0;tautau_est_lby2_av=0.0;tautau_est_lby3_av=0.0;tautau_est_lby6_av=0.0;psistarpsi_et_lby2_av=0.0;psistarpsi_et_lby3_av=0.0;psistarpsi_et_lby6_av=0.0;psistarpsi_stat_lby2_av=0.0;psistarpsi_stat_lby3_av=0.0;psistarpsi_stat_lby6_av=0.0;psistarpsi_est_lby2_av=0.0;psistarpsi_est_lby3_av=0.0;psistarpsi_est_lby6_av=0.0;psi2starpsi2_est_lby2_av=0.0;psi2starpsi2_est_lby3_av=0.0;psi2starpsi2_est_lby6_av=0.0;psi3starpsi3_est_lby2_av=0.0;psi3starpsi3_est_lby3_av=0.0;psi3starpsi3_est_lby6_av=0.0;
psi3tau_est_lby2_av=psi3tau_est_lby3_av=psi3tau_est_lby6_av=0.0;
impsi_av=repsi_av=0.0;
  n_measure=0;
  completed_bins=0;
for (i=0;i<63;i++) {
hist_phase[i]=0; }
  return;
}
