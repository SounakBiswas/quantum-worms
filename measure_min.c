#include"global.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define  next_nbr(ii,jj,xx,yy)  ((xx+ii+lx)%lx + ((yy+jj+ly)%ly)*lx)

void  bin(  double, double, double, double, double, double); 


void measure()
{
  double energy_conf;
 // double fac=beta/((double)(opstr_l*(opstr_l+1)));
 // double fac1=1.0/((double)opstr_l);
 double fac=beta/((double)(n_niop*(n_niop+1)));
 double fac1=beta/((double)(n_niop));
  int p_counter; //index of current niop
  int last_flip;
  int n_flip[NSITES];
  int last_nio;
double sigmaz_q_stat_conf;
double sigmaz_q_et_conf;
double sigmaz_0_et_conf;
double sigmaz_0_stat_conf;
double sigmax_q_stat_conf;
/* Variables/Buffers to be sent to bin*****/
 


/***Other Buffers,for static guys*******/
double resigmaz_q_sumrp,imsigmaz_q_sumrp,sigmaz_0_sumrp;
double resigmax_q_stat,imsigmax_q_stat;
/*************************/
double resigmaz_q_p,imsigmaz_q_p,sigmaz_0_p;// summed over r, at a particular time slice
int i,j,k,p; 
p_counter=last_flip=0;
/*Initialization of Buffers*****/
n_measure++;
sigmaz_q_stat_conf=sigmaz_q_et_conf=sigmaz_0_et_conf=sigmaz_0_stat_conf=sigmax_q_stat_conf=0;
resigmaz_q_sumrp=imsigmaz_q_sumrp=sigmaz_0_sumrp=0;
 energy_conf=-(double)n_niop/beta;
/******/
/*Initialization of Feeders*/
  resigmaz_q_p=imsigmaz_q_p=sigmaz_0_p=0;
  for (i=0;i<nsites;i++){
    resigmaz_q_p+=sigma[i]*Cos[sublat[i]];
    imsigmaz_q_p+=sigma[i]*Sin[sublat[i]];
    sigmaz_0_p+=sigma[i];
    last_flip=n_flip[i]=0;
  }

/**Pass through the operator  string to add to buffers*/
 double change;
 int site1;
  n_flipop=0;
 for(p=0;p<opstr_l;p++){
   if(opstr[p]!=-1){
   p_counter++;
    last_nio=p; }
  if(opstr[p]>=ntrianglesandsites){
   sigmaz_q_et_conf+=(p_counter-last_flip)*(resigmaz_q_p*resigmaz_q_p+imsigmaz_q_p*imsigmaz_q_p);
   sigmaz_0_et_conf+=(p_counter-last_flip)*(sigmaz_0_p*sigmaz_0_p);
   resigmaz_q_sumrp+=(p_counter-last_flip)*(resigmaz_q_p);   
   imsigmaz_q_sumrp+=(p_counter-last_flip)*imsigmaz_q_p;
   sigmaz_0_sumrp+=(p_counter-last_flip)*sigmaz_0_p; 
   /*******************/ 
   last_flip=p_counter;
   n_flipop++;
   site1=opstr[p]-ntrianglesandsites;
   n_flip[site1]++;
/* Editing the feeders*/
   
 //feeders for global qtys
   resigmaz_q_p-=(2*sigma[site1]*Cos[sublat[site1]]);
   imsigmaz_q_p-=(2*sigma[site1]*Sin[sublat[site1]]); 
   sigmaz_0_p-=(2*sigma[site1]);
//propagate the spin state
   sigma[site1]*=-1;


  }

 }
/***Addition to the buffer at the end of opstring*******/
 if(opstr[last_nio]<ntrianglesandsites){
  
   sigmaz_q_et_conf+=(p_counter-last_flip)*(resigmaz_q_p*resigmaz_q_p+imsigmaz_q_p*imsigmaz_q_p);
   sigmaz_0_et_conf+=(p_counter-last_flip)*(sigmaz_0_p*sigmaz_0_p);
   resigmaz_q_sumrp+=(p_counter-last_flip)*(resigmaz_q_p);   
   imsigmaz_q_sumrp+=(p_counter-last_flip)*imsigmaz_q_p;
   sigmaz_0_sumrp+=(p_counter-last_flip)*sigmaz_0_p; 

 }
 /*Final loop over sites to get the equaltime operators,the fake quantities, and also do the work for sigmax correlator*/

 resigmax_q_stat=imsigmax_q_stat=0;
 for(i=0;i<nsites;i++) {
   site1=i;
   resigmax_q_stat+=n_flip[site1]*Cos[sublat[site1]];
   imsigmax_q_stat+=n_flip[site1]*Sin[sublat[site1]];
 }  
 
 sigmaz_q_stat_conf=sigmaz_q_et_conf+(resigmaz_q_sumrp*resigmaz_q_sumrp+imsigmaz_q_sumrp*imsigmaz_q_sumrp);
 sigmaz_0_stat_conf=sigmaz_0_et_conf+(sigmaz_0_sumrp*sigmaz_0_sumrp);
 sigmax_q_stat_conf=resigmax_q_stat*resigmax_q_stat+imsigmax_q_stat*imsigmax_q_stat-n_flipop;
  /// Normalization
  sigmax_q_stat_conf = sigmax_q_stat_conf*(1.0/(beta*h*h));
  sigmaz_q_stat_conf = sigmaz_q_stat_conf*fac;
  sigmaz_q_et_conf = sigmaz_q_et_conf*fac1;
  sigmaz_0_stat_conf =sigmaz_0_stat_conf*fac;
  sigmaz_0_et_conf = sigmaz_0_et_conf*fac1;
  
 //auto_corr(modsigmax_q_int,modsigmaz_q_int,modsigmaz_0_int,psistarpsi_et_lby6,psistarpsi_stat_lby6);

  bin(energy_conf,sigmax_q_stat_conf,sigmaz_q_stat_conf,sigmaz_q_et_conf,sigmaz_0_stat_conf,sigmaz_0_et_conf);


  return;


}


