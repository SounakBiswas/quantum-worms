#include"global.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define  next_nbr(ii,jj,xx,yy)  ((xx+ii+lx)%lx + ((yy+jj+ly)%ly)*lx)

void  bin(  double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );

void auto_corr(double op1,double op2,double op3,double op4,double op5);


void measure()
{
  double energy_conf;
  double fac=beta/((double)(n_niop*(n_niop+1)));
  double fac1=1.0/((double)n_niop);
  //double fac=beta/((double)(n_niop*(n_niop+1)));
  //double fac1=1.0/((double)n_niop);
  
  int p_counter; //index of current niop
  int last_flip;
  int last_change_psi[NSITES];
  int n_flip[NSITES];
double sigmaz_q_stat_conf;
double sigmaz_q_et_conf;
double sigmaz_0_et_conf;
double sigmaz_0_stat_conf;
double sigmax_q_stat_conf;
/* Variables/Buffers to be sent to bin*****/
//Tau Honest
double tautau_et_lby2_conf,tautau_et_lby3_conf,tautau_et_lby6_conf;
double tautau_stat_lby2_conf,tautau_stat_lby3_conf,tautau_stat_lby6_conf;
double tautau_est_lby2_conf,tautau_est_lby3_conf,tautau_est_lby6_conf;
 
//Psi honest
double psistarpsi_et_lby2_conf,psistarpsi_et_lby3_conf,psistarpsi_et_lby6_conf;
double psistarpsi_stat_lby2_conf,psistarpsi_stat_lby3_conf,psistarpsi_stat_lby6_conf;
double psistarpsi_est_lby2_conf,psistarpsi_est_lby3_conf,psistarpsi_est_lby6_conf;
//Fake guys
double psi2starpsi2_est_lby2_conf,psi2starpsi2_est_lby3_conf,psi2starpsi2_est_lby6_conf;
double psi3starpsi3_est_lby2_conf,psi3starpsi3_est_lby3_conf,psi3starpsi3_est_lby6_conf;
double psi3tau_est_lby2_conf,psi3tau_est_lby3_conf,psi3tau_est_lby6_conf;

double phase_psi_conf;
double repsi_sumrp_conf,impsi_sumrp_conf;
/***************************************************/

/***Other Buffers,for static guys*******/
double repsi_sump[NSITES],impsi_sump[NSITES],tau_sump[NSITES];double resigmaz_q_sumrp,imsigmaz_q_sumrp,sigmaz_0_sumrp;
double resigmax_q_stat,imsigmax_q_stat;
/*************************/

/**current values of things, ie at a particular time slice, feeders to buffers***/
double repsi_p[NSITES],impsi_p[NSITES],tau_p[NSITES];
double psistarpsi_lby2_p,psistarpsi_lby3_p,psistarpsi_lby6_p,tautau_lby2_p,tautau_lby3_p,tautau_lby6_p;
/*********************/
double resigmaz_q_p,imsigmaz_q_p,sigmaz_0_p;// summed over r, at a particular time slice
int i,j,k,p; 
p_counter=last_flip=0;
/*Initialization of Buffers*****/
n_measure++;
sigmaz_q_stat_conf=sigmaz_q_et_conf=sigmaz_0_et_conf=sigmaz_0_stat_conf=sigmax_q_stat_conf=tautau_et_lby2_conf=tautau_et_lby3_conf=tautau_et_lby6_conf=psistarpsi_et_lby2_conf=psistarpsi_et_lby3_conf=psistarpsi_et_lby6_conf=0;
resigmaz_q_sumrp=imsigmaz_q_sumrp=sigmaz_0_sumrp=0;
repsi_sumrp_conf=impsi_sumrp_conf=0;
for(i=0;i<nsites;i++)
 repsi_sump[i]=impsi_sump[i]=tau_sump[i]=0;

 phase_psi_conf=0;
 energy_conf=-(double)n_niop/beta;
/******/
/*Initialization of Feeders*/
  int up_trian;
  int sitelby2,sitelby3,sitelby6,siteminuslby2,siteminuslby3,siteminuslby6;
  int sitea,siteb,sitec,site;
  resigmaz_q_p=imsigmaz_q_p=sigmaz_0_p=psistarpsi_lby2_p=psistarpsi_lby3_p=psistarpsi_lby6_p=tautau_lby2_p=tautau_lby3_p=tautau_lby6_p=0;
  for (i=0;i<nsites;i++){
    up_trian=2*i;
    sitea=triangle[up_trian][0];siteb=triangle[up_trian][1];sitec=triangle[up_trian][2];
    repsi_p[i]=sigma[sitea]*Cos[0]+sigma[siteb]*Cos[1]+sigma[sitec]*Cos[2];
    impsi_p[i]=sigma[sitea]*Sin[0]+sigma[siteb]*Sin[1]+sigma[sitec]*Sin[2];   
    tau_p[i]=sigma[sitea]+sigma[siteb]+sigma[sitec];
    resigmaz_q_p+=sigma[i]*Cos[sublat[i]];
    imsigmaz_q_p+=sigma[i]*Sin[sublat[i]];
    sigmaz_0_p+=sigma[i];
    last_flip=last_change_psi[i]=n_flip[i]=0;
  }
  for(i=0;i<nsites;i++){
    
    j=i/lx; k=i%lx; 
    sitelby2=next_nbr(k,j,lx/2,0);
    sitelby3=next_nbr(k,j,lx/3,0);     
    sitelby6=next_nbr(k,j,lx/6,0);
    psistarpsi_lby2_p+=repsi_p[i]*repsi_p[sitelby2]+impsi_p[i]*impsi_p[sitelby2];
    psistarpsi_lby3_p+=repsi_p[i]*repsi_p[sitelby3]+impsi_p[i]*impsi_p[sitelby3];
    psistarpsi_lby6_p+=repsi_p[i]*repsi_p[sitelby6]+impsi_p[i]*impsi_p[sitelby6];
    tautau_lby2_p+=tau_p[i]*tau_p[sitelby2];
    tautau_lby3_p+=tau_p[i]*tau_p[sitelby3];
    tautau_lby6_p+=tau_p[i]*tau_p[sitelby6];
  }

 /********************************************/
/**Pass through the operator  string to add to buffers*/
 int site1,site2,site3;//other sites affected when spin filps at site1
int site1x,site1y,site2x,site2y,site3x,site3y;
 double change;
 double test =0;
  n_flipop=0;
 for(p=0;p<opstr_l;p++){
 if(opstr[p]!=-1)
   p_counter++;
   
  if(opstr[p]>=ntrianglesandsites){
   site1=opstr[p]-ntrianglesandsites;
   site2=nbr[site1][3];
   site3=nbr[site1][4]; 
   site1x=site1%lx;site1y=site1/lx;
   site2x=site2%lx;site2y=site2/lx;
   site3x=site3%lx;site3y=site3/lx;
//   printf("%d\t%d\t%d\t%d\t%d\n",site1,next_nbr(site1x,site1y,lx/2,0),next_nbr(site1x,site1y,-lx/2,0),next_nbr(site1x,site1y,lx/6,0),next_nbr(site1x,site1y,-lx/6,0));
  // getchar();
   /***add to buffers**/
   psistarpsi_et_lby2_conf+=(p_counter-last_flip)*psistarpsi_lby2_p;
   psistarpsi_et_lby3_conf+=(p_counter-last_flip)*psistarpsi_lby3_p;
   psistarpsi_et_lby6_conf+=(p_counter-last_flip)*psistarpsi_lby6_p;
   tautau_et_lby2_conf+=(p_counter-last_flip)*tautau_lby2_p;
   tautau_et_lby3_conf+=(p_counter-last_flip)*tautau_lby3_p;
   tautau_et_lby6_conf+=(p_counter-last_flip)*tautau_lby6_p;
   sigmaz_q_et_conf+=(p_counter-last_flip)*(resigmaz_q_p*resigmaz_q_p+imsigmaz_q_p*imsigmaz_q_p);
   sigmaz_0_et_conf+=(p_counter-last_flip)*(sigmaz_0_p*sigmaz_0_p);

   repsi_sump[site1]+=(p_counter-last_change_psi[site1])*repsi_p[site1];  // changing psi and tau for three sites
   impsi_sump[site1]+=(p_counter-last_change_psi[site1])*impsi_p[site1];
   tau_sump[site1]+=(p_counter-last_change_psi[site1])*tau_p[site1];
   repsi_sump[site2]+=(p_counter-last_change_psi[site2])*repsi_p[site2];
   impsi_sump[site2]+=(p_counter-last_change_psi[site2])*impsi_p[site2];
   tau_sump[site2]+=(p_counter-last_change_psi[site2])*tau_p[site2];
   repsi_sump[site3]+=(p_counter-last_change_psi[site3])*repsi_p[site3];
   impsi_sump[site3]+=(p_counter-last_change_psi[site3])*impsi_p[site3];
   tau_sump[site3]+=(p_counter-last_change_psi[site3])*tau_p[site3];

   resigmaz_q_sumrp+=(p_counter-last_flip)*(resigmaz_q_p);   
   imsigmaz_q_sumrp+=(p_counter-last_flip)*imsigmaz_q_p;
   sigmaz_0_sumrp+=(p_counter-last_flip)*sigmaz_0_p; 
   /*******************/ 
   last_change_psi[site1]=last_change_psi[site2]=last_change_psi[site3]=last_flip=p_counter;
   n_flip[site1]++;
   n_flipop++;
/* Editing the feeders*/
   
   psistarpsi_lby2_p-=(repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,lx/2,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,lx/2,0)]+repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,-lx/2,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,-lx/2,0)]);
   psistarpsi_lby2_p-=(repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,lx/2,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,lx/2,0)]+repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,-lx/2,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,-lx/2,0)]);
   psistarpsi_lby2_p-=(repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,lx/2,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,lx/2,0)]+repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,-lx/2,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,-lx/2,0)]);
   
   psistarpsi_lby3_p-=(repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,lx/3,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,lx/3,0)]+repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,-lx/3,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,-lx/3,0)]);
   psistarpsi_lby3_p-=(repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,lx/3,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,lx/3,0)]+repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,-lx/3,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,-lx/3,0)]);
   psistarpsi_lby3_p-=(repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,lx/3,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,lx/3,0)]+repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,-lx/3,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,-lx/3,0)]);
    
   psistarpsi_lby6_p-=(repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,lx/6,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,lx/6,0)]+repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,-lx/6,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,-lx/6,0)]);
   psistarpsi_lby6_p-=(repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,lx/6,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,lx/6,0)]+repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,-lx/6,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,-lx/6,0)]);
   psistarpsi_lby6_p-=(repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,lx/6,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,lx/6,0)]+repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,-lx/6,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,-lx/6,0)]);
   
   tautau_lby2_p-=(tau_p[site1]*tau_p[next_nbr(site1x,site1y,lx/2,0)]+tau_p[site1]*tau_p[next_nbr(site1x,site1y,-lx/2,0)]);
   tautau_lby2_p-=(tau_p[site2]*tau_p[next_nbr(site2x,site2y,lx/2,0)]+tau_p[site2]*tau_p[next_nbr(site2x,site2y,-lx/2,0)]);
   tautau_lby2_p-=(tau_p[site3]*tau_p[next_nbr(site3x,site3y,lx/2,0)]+tau_p[site3]*tau_p[next_nbr(site3x,site3y,-lx/2,0)]);
  
   tautau_lby3_p-=(tau_p[site1]*tau_p[next_nbr(site1x,site1y,lx/3,0)]+tau_p[site1]*tau_p[next_nbr(site1x,site1y,-lx/3,0)]);
   tautau_lby3_p-=(tau_p[site2]*tau_p[next_nbr(site2x,site2y,lx/3,0)]+tau_p[site2]*tau_p[next_nbr(site2x,site2y,-lx/3,0)]);
   tautau_lby3_p-=(tau_p[site3]*tau_p[next_nbr(site3x,site3y,lx/3,0)]+tau_p[site3]*tau_p[next_nbr(site3x,site3y,-lx/3,0)]);
   
   tautau_lby6_p-=(tau_p[site1]*tau_p[next_nbr(site1x,site1y,lx/6,0)]+tau_p[site1]*tau_p[next_nbr(site1x,site1y,-lx/6,0)]);
   tautau_lby6_p-=(tau_p[site2]*tau_p[next_nbr(site2x,site2y,lx/6,0)]+tau_p[site2]*tau_p[next_nbr(site2x,site2y,-lx/6,0)]);
   tautau_lby6_p-=(tau_p[site3]*tau_p[next_nbr(site3x,site3y,lx/6,0)]+tau_p[site3]*tau_p[next_nbr(site3x,site3y,-lx/6,0)]);
 //***********************
 //feeders for global qtys
   resigmaz_q_p-=2*sigma[site1]*Cos[sublat[site1]];
   imsigmaz_q_p-=2*sigma[site1]*Sin[sublat[site1]]; 
   sigmaz_0_p-=2*sigma[site1];
//feeders for static qtys
   change=-2*Cos[sublat[site1]]*sigma[site1];  
   repsi_p[site1]+=change;repsi_p[site2]+=change;repsi_p[site3]+=change;
   change=-2*sigma[site1]*Sin[sublat[site1]];
   impsi_p[site1]+=change;impsi_p[site2]+=change;impsi_p[site3]+=change;
   change=-2*sigma[site1];
   tau_p[site1]+=change;tau_p[site2]+=change;tau_p[site3]+=change;

//propagate the spin state
   sigma[site1]*=-1;


//finally change the feeders for the equal time correlators
   psistarpsi_lby2_p+=(repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,lx/2,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,lx/2,0)]+repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,-lx/2,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,-lx/2,0)]);
   psistarpsi_lby2_p+=(repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,lx/2,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,lx/2,0)]+repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,-lx/2,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,-lx/2,0)]);
   psistarpsi_lby2_p+=(repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,lx/2,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,lx/2,0)]+repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,-lx/2,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,-lx/2,0)]);
   
   psistarpsi_lby3_p+=(repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,lx/3,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,lx/3,0)]+repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,-lx/3,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,-lx/3,0)]);
   psistarpsi_lby3_p+=(repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,lx/3,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,lx/3,0)]+repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,-lx/3,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,-lx/3,0)]);
   psistarpsi_lby3_p+=(repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,lx/3,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,lx/3,0)]+repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,-lx/3,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,-lx/3,0)]);
    
   psistarpsi_lby6_p+=(repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,lx/6,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,lx/6,0)]+repsi_p[site1]*repsi_p[next_nbr(site1x,site1y,-lx/6,0)]+impsi_p[site1]*impsi_p[next_nbr(site1x,site1y,-lx/6,0)]);
   psistarpsi_lby6_p+=(repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,lx/6,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,lx/6,0)]+repsi_p[site2]*repsi_p[next_nbr(site2x,site2y,-lx/6,0)]+impsi_p[site2]*impsi_p[next_nbr(site2x,site2y,-lx/6,0)]);
   psistarpsi_lby6_p+=(repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,lx/6,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,lx/6,0)]+repsi_p[site3]*repsi_p[next_nbr(site3x,site3y,-lx/6,0)]+impsi_p[site3]*impsi_p[next_nbr(site3x,site3y,-lx/6,0)]);
   
   tautau_lby2_p+=(tau_p[site1]*tau_p[next_nbr(site1x,site1y,lx/2,0)]+tau_p[site1]*tau_p[next_nbr(site1x,site1y,-lx/2,0)]);
   tautau_lby2_p+=(tau_p[site2]*tau_p[next_nbr(site2x,site2y,lx/2,0)]+tau_p[site2]*tau_p[next_nbr(site2x,site2y,-lx/2,0)]);
   tautau_lby2_p+=(tau_p[site3]*tau_p[next_nbr(site3x,site3y,lx/2,0)]+tau_p[site3]*tau_p[next_nbr(site3x,site3y,-lx/2,0)]);
  
   tautau_lby3_p+=(tau_p[site1]*tau_p[next_nbr(site1x,site1y,lx/3,0)]+tau_p[site1]*tau_p[next_nbr(site1x,site1y,-lx/3,0)]);
   tautau_lby3_p+=(tau_p[site2]*tau_p[next_nbr(site2x,site2y,lx/3,0)]+tau_p[site2]*tau_p[next_nbr(site2x,site2y,-lx/3,0)]);
   tautau_lby3_p+=(tau_p[site3]*tau_p[next_nbr(site3x,site3y,lx/3,0)]+tau_p[site3]*tau_p[next_nbr(site3x,site3y,-lx/3,0)]);
   
   tautau_lby6_p+=(tau_p[site1]*tau_p[next_nbr(site1x,site1y,lx/6,0)]+tau_p[site1]*tau_p[next_nbr(site1x,site1y,-lx/6,0)]);
   tautau_lby6_p+=(tau_p[site2]*tau_p[next_nbr(site2x,site2y,lx/6,0)]+tau_p[site2]*tau_p[next_nbr(site2x,site2y,-lx/6,0)]);
   tautau_lby6_p+=(tau_p[site3]*tau_p[next_nbr(site3x,site3y,lx/6,0)]+tau_p[site3]*tau_p[next_nbr(site3x,site3y,-lx/6,0)]);
  }

 }
/***Addition to the buffer at the end of opstring*******/
 
  
   psistarpsi_et_lby2_conf+=(p_counter-last_flip)*psistarpsi_lby2_p;
   psistarpsi_et_lby3_conf+=(p_counter-last_flip)*psistarpsi_lby3_p;
   psistarpsi_et_lby6_conf+=(p_counter-last_flip)*psistarpsi_lby6_p;
   tautau_et_lby2_conf+=(p_counter-last_flip)*tautau_lby2_p;
   tautau_et_lby3_conf+=(p_counter-last_flip)*tautau_lby3_p;
   tautau_et_lby6_conf+=(p_counter-last_flip)*tautau_lby6_p;
   sigmaz_q_et_conf+=(p_counter-last_flip)*(resigmaz_q_p*resigmaz_q_p+imsigmaz_q_p*imsigmaz_q_p);
   sigmaz_0_et_conf+=(p_counter-last_flip)*(sigmaz_0_p*sigmaz_0_p);
   resigmaz_q_sumrp+=(p_counter-last_flip)*(resigmaz_q_p);   
   imsigmaz_q_sumrp+=(p_counter-last_flip)*imsigmaz_q_p;
   sigmaz_0_sumrp+=(p_counter-last_flip)*sigmaz_0_p; 
   for(i=0;i<nsites;i++){
    repsi_sump[i]+=(p_counter-last_change_psi[i])*repsi_p[i];  // changing psi and tau for three sites
    impsi_sump[i]+=(p_counter-last_change_psi[i])*impsi_p[i];
    tau_sump[i]+=(p_counter-last_change_psi[i])*tau_p[i];
   }

 
 /*Final loop over sites to get the static operators,the fake quantities, and also do the work for sigmax correlator*/

 psistarpsi_stat_lby2_conf=psistarpsi_stat_lby3_conf=psistarpsi_stat_lby6_conf=tautau_stat_lby2_conf=tautau_stat_lby3_conf=tautau_stat_lby6_conf=0.0;
psi2starpsi2_est_lby2_conf=psi2starpsi2_est_lby3_conf=psi2starpsi2_est_lby6_conf=psi3starpsi3_est_lby2_conf=psi3starpsi3_est_lby3_conf=psi3starpsi3_est_lby6_conf= psi3tau_est_lby2_conf=psi3tau_est_lby3_conf=psi3tau_est_lby6_conf=tautau_est_lby2_conf=tautau_est_lby3_conf=tautau_est_lby6_conf= psistarpsi_est_lby2_conf=psistarpsi_est_lby3_conf=psistarpsi_est_lby6_conf=0;
 resigmax_q_stat=imsigmax_q_stat=0;
 repsi_sumrp_conf=impsi_sumrp_conf=0;
 double a,b,c,d,a2,b2,c2,d2,c3,d3,a3,b3; //Psi at p,summed over r,is a +ib
 double e,f,e2,f2,e3,f3;
 for(i=0;i<nsites;i++) {
   site1=i;
   site1x=site1%lx;site1y=site1/lx;
   sitelby2=next_nbr(site1x,site1y,lx/2,0);
   sitelby3=next_nbr(site1x,site1y,lx/3,0);
   sitelby6=next_nbr(site1x,site1y,lx/6,0);
   
   a=repsi_sump[site1];
   b=impsi_sump[site1];
   e=tau_sump[site1];
   c=repsi_sump[sitelby2];
   d=impsi_sump[sitelby2];
   f=tau_sump[sitelby2];
   a2=a*a;a3=a2*a;b2=b*b;b3=b2*b;
   c2=c*c;c3=c2*c;d2=d*d;d3=d2*d;
   psistarpsi_est_lby2_conf+=(a*c+b*d);
   tautau_est_lby2_conf+=e*f;
   psi2starpsi2_est_lby2_conf+=(a2-b2)*(c2-d2)+4*a*b*c*d;
   psi3starpsi3_est_lby2_conf+=(a3-3*a*b2)*(c3-3*c*d2)+(3*a2*b-b3)*(3*c2*d-d3);
   psi3tau_est_lby2_conf+=(a3-3*a*b2)*f+e*(c3-3*c*d2);
  
   c=repsi_sump[sitelby3];
   d=impsi_sump[sitelby3]; 
   f=tau_sump[sitelby3];
   c2=c*c;c3=c2*c;d2=d*d;d3=d2*d;
   psistarpsi_est_lby3_conf+=(a*c+b*d);
   tautau_est_lby3_conf+=e*f;
   psi2starpsi2_est_lby3_conf+=(a2-b2)*(c2-d2)+4*a*b*c*d;
   psi3starpsi3_est_lby3_conf+=(a3-3*a*b2)*(c3-3*c*d2)+(3*a2*b-b3)*(3*c2*d-d3);
   psi3tau_est_lby3_conf+=(a3-3*a*b2)*f+e*(c3-3*c*d2);
   
   c=repsi_sump[sitelby6];
   d=impsi_sump[sitelby6];
   f=tau_sump[sitelby6];
   c2=c*c;c3=c2*c;d2=d*d;d3=d2*d;
   psistarpsi_est_lby6_conf+=(a*c+b*d);
   tautau_est_lby6_conf+=e*f;
   psi2starpsi2_est_lby6_conf+=(a2-b2)*(c2-d2)+4*a*b*c*d;
   psi3starpsi3_est_lby6_conf+=(a3-3*a*b2)*(c3-3*c*d2)+(3*a2*b-b3)*(3*c2*d-d3);
   psi3tau_est_lby6_conf+=(a3-3*a*b2)*f+e*(c3-3*c*d2);
   
   
   resigmax_q_stat+=n_flip[site1]*Cos[sublat[site1]];
   imsigmax_q_stat+=n_flip[site1]*Sin[sublat[site1]];
   repsi_sumrp_conf+=repsi_sump[site1];
   impsi_sumrp_conf+=impsi_sump[site1];
 }  
 
 psistarpsi_stat_lby2_conf=psistarpsi_et_lby2_conf+psistarpsi_est_lby2_conf;
 tautau_stat_lby2_conf=tautau_et_lby2_conf+tautau_est_lby2_conf;
 psistarpsi_stat_lby3_conf=psistarpsi_et_lby3_conf+psistarpsi_est_lby3_conf;
 tautau_stat_lby3_conf=tautau_et_lby3_conf+tautau_est_lby3_conf;
 psistarpsi_stat_lby6_conf=psistarpsi_et_lby6_conf+psistarpsi_est_lby6_conf;
 tautau_stat_lby6_conf=tautau_et_lby6_conf+tautau_est_lby6_conf;
 sigmaz_q_stat_conf=sigmaz_q_et_conf+(resigmaz_q_sumrp*resigmaz_q_sumrp+imsigmaz_q_sumrp*imsigmaz_q_sumrp);
 sigmaz_0_stat_conf=sigmaz_0_et_conf+(sigmaz_0_sumrp*sigmaz_0_sumrp);
 sigmax_q_stat_conf=resigmax_q_stat*resigmax_q_stat+imsigmax_q_stat*imsigmax_q_stat-n_flipop;
 if (fabs(repsi_sumrp_conf) > 0.000001) {
    phase_psi_conf = atan(impsi_sumrp_conf/repsi_sumrp_conf);  
  }
  else {
    phase_psi_conf = copysign(piby2,impsi_sumrp_conf) * copysign(1.0,repsi_sumrp_conf);
  }
  // Histogram of phase of psi 
   hist_phase[(int)((phase_psi_conf+pi)*10)]++;
  // histogram is being printed in bin.c

  double repsi_conf=repsi_sumrp_conf*fac1;
  double impsi_conf=impsi_sumrp_conf*fac1;

  /// Normalization
  sigmax_q_stat_conf = sigmax_q_stat_conf*(1.0/(beta*h*h));
  sigmaz_q_stat_conf = sigmaz_q_stat_conf*fac;
  sigmaz_q_et_conf = sigmaz_q_et_conf*fac1;
  sigmaz_0_stat_conf =sigmaz_0_stat_conf*fac;
  sigmaz_0_et_conf = sigmaz_0_et_conf*fac1;
  
  psistarpsi_et_lby2_conf *= fac1;
  psistarpsi_et_lby3_conf *= fac1;
  psistarpsi_et_lby6_conf *= fac1;
  psistarpsi_stat_lby2_conf *= fac;
  psistarpsi_stat_lby3_conf *= fac;
  psistarpsi_stat_lby6_conf *= fac;
  
  tautau_et_lby2_conf *= fac1;
  tautau_et_lby3_conf *= fac1;
  tautau_et_lby6_conf *= fac1;
  tautau_stat_lby2_conf *= fac;
  tautau_stat_lby3_conf *= fac;
  tautau_stat_lby6_conf *= fac;
  
  psistarpsi_est_lby2_conf*=fac1*fac1;
  psistarpsi_est_lby3_conf*=fac1*fac1;
  psistarpsi_est_lby6_conf*=fac1*fac1;
  psi2starpsi2_est_lby2_conf*=fac1*fac1*fac1*fac1;
  psi2starpsi2_est_lby3_conf*=fac1*fac1*fac1*fac1;;
  psi2starpsi2_est_lby6_conf*=fac1*fac1*fac1*fac1;
  psi3starpsi3_est_lby2_conf*=fac1*fac1*fac1*fac1*fac1*fac1;
  psi3starpsi3_est_lby3_conf*=fac1*fac1*fac1*fac1*fac1*fac1;
  psi3starpsi3_est_lby6_conf*=fac1*fac1*fac1*fac1*fac1*fac1;
  psi3tau_est_lby2_conf*=fac1*fac1*fac1*fac1;
  psi3tau_est_lby3_conf*=fac1*fac1*fac1*fac1;
  psi3tau_est_lby6_conf*=fac1*fac1*fac1*fac1;
  tautau_est_lby2_conf*=fac1*fac1;
  tautau_est_lby3_conf*=fac1*fac1;
  tautau_est_lby6_conf*=fac1*fac1;

 auto_corr(sigmax_q_stat_conf,sigmaz_q_stat_conf,sigmaz_0_stat_conf,psistarpsi_et_lby6_conf,psistarpsi_stat_lby6_conf);

  bin(energy_conf,sigmax_q_stat_conf,sigmaz_q_stat_conf,sigmaz_q_et_conf,sigmaz_0_stat_conf,sigmaz_0_et_conf,psistarpsi_et_lby2_conf,psistarpsi_et_lby3_conf,psistarpsi_et_lby6_conf,psistarpsi_stat_lby2_conf,psistarpsi_stat_lby3_conf,psistarpsi_stat_lby6_conf,
tautau_et_lby2_conf,tautau_et_lby3_conf,tautau_et_lby6_conf,tautau_stat_lby2_conf,tautau_stat_lby3_conf,tautau_stat_lby6_conf,
psistarpsi_est_lby2_conf,psistarpsi_est_lby3_conf,psistarpsi_est_lby6_conf,psi2starpsi2_est_lby2_conf,psi2starpsi2_est_lby3_conf,psi2starpsi2_est_lby6_conf,psi3starpsi3_est_lby2_conf,psi3starpsi3_est_lby3_conf,psi3starpsi3_est_lby6_conf,
tautau_est_lby2_conf,tautau_est_lby3_conf,tautau_est_lby6_conf,
psi3tau_est_lby2_conf,psi3tau_est_lby3_conf,psi3tau_est_lby6_conf,phase_psi_conf,repsi_conf,impsi_conf);


  return;


}


