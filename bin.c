#include<stdio.h>
#include "global.h"
void bin(double energy,double sigmax_q_stat,double sigmaz_q_stat,double sigmaz_q_et,double sigmaz_0_stat,double sigmaz_0_et,double psistarpsi_et_lby2,double psistarpsi_et_lby3,double psistarpsi_et_lby6,double psistarpsi_stat_lby2,double psistarpsi_stat_lby3,double psistarpsi_stat_lby6,
    double tautau_et_lby2,double tautau_et_lby3,double tautau_et_lby6,double tautau_stat_lby2,double tautau_stat_lby3,double tautau_stat_lby6,
    double psistarpsi_est_lby2,double psistarpsi_est_lby3,double psistarpsi_est_lby6,double psi2starpsi2_est_lby2,double psi2starpsi2_est_lby3,double psi2starpsi2_est_lby6,double psi3starpsi3_est_lby2,double psi3starpsi3_est_lby3,double psi3starpsi3_est_lby6,
    double tautau_est_lby2,double tautau_est_lby3,double tautau_est_lby6,
    double psi3tau_est_lby2,double psi3tau_est_lby3,double psi3tau_est_lby6,double phase_psi,double repsi,double impsi){

int i;

  FILE *binf, *cumf,*histf, *clustf;

  static long num_call_bin=0;

  static double  energy_bin;
  static double sigmax_q_stat_bin;
  static double sigmaz_q_stat_bin;
  static double sigmaz_q_et_bin;
  static double sigmaz_0_stat_bin;
  static double sigmaz_0_et_bin;


  static double psistarpsi_et_lby2_bin,psistarpsi_et_lby3_bin,psistarpsi_et_lby6_bin;
  static double psistarpsi_stat_lby2_bin,psistarpsi_stat_lby3_bin,psistarpsi_stat_lby6_bin;
  static double tautau_et_lby2_bin,tautau_et_lby3_bin,tautau_et_lby6_bin;
  static double tautau_stat_lby2_bin,tautau_stat_lby3_bin,tautau_stat_lby6_bin;
  static double psistarpsi_est_lby2_bin,psistarpsi_est_lby3_bin,psistarpsi_est_lby6_bin;
  static double psi2starpsi2_est_lby2_bin,psi2starpsi2_est_lby3_bin,psi2starpsi2_est_lby6_bin;
  static double psi3starpsi3_est_lby2_bin,psi3starpsi3_est_lby3_bin,psi3starpsi3_est_lby6_bin;
  static double tautau_est_lby2_bin,tautau_est_lby3_bin,tautau_est_lby6_bin;
  static double psi3tau_est_lby2_bin,psi3tau_est_lby3_bin,psi3tau_est_lby6_bin;
  static double phase_psi_bin;
  static double repsi_bin;
  static double impsi_bin;


  //initialise bin variables

  if(!(num_call_bin % binsize))
  {
    energy_bin=0.0;
    sigmax_q_stat_bin=0.0;
    sigmaz_q_stat_bin=0.0;
    sigmaz_q_et_bin=0.0;
    sigmaz_0_stat_bin=0.0;
    sigmaz_0_et_bin=0.0;


    psistarpsi_et_lby2_bin=psistarpsi_et_lby3_bin=psistarpsi_et_lby6_bin=0.0;
    psistarpsi_stat_lby2_bin=psistarpsi_stat_lby3_bin=psistarpsi_stat_lby6_bin=0.0;
    tautau_et_lby2_bin=tautau_et_lby3_bin=tautau_et_lby6_bin=0.0;
    tautau_stat_lby2_bin=tautau_stat_lby3_bin=tautau_stat_lby6_bin=0.0;
    psistarpsi_est_lby2_bin=psistarpsi_est_lby3_bin=psistarpsi_est_lby6_bin=0.0;
    psi2starpsi2_est_lby2_bin=psi2starpsi2_est_lby3_bin=psi2starpsi2_est_lby6_bin=0.0;
    psi3starpsi3_est_lby2_bin=psi3starpsi3_est_lby3_bin=psi3starpsi3_est_lby6_bin=0.0;
    tautau_est_lby2_bin=tautau_est_lby3_bin=tautau_est_lby6_bin=0.0;
    psi3tau_est_lby2_bin=psi3tau_est_lby3_bin=psi3tau_est_lby6_bin=0.0;
    phase_psi_bin=0.0;
    repsi_bin=0.0;
    impsi_bin=0.0;

  }

  num_call_bin++;	

  //incrementing bin variables

  energy_bin+=energy;
  sigmax_q_stat_bin+=sigmax_q_stat;
  sigmaz_q_stat_bin+= sigmaz_q_stat;
  sigmaz_q_et_bin+=sigmaz_q_et;
  sigmaz_0_stat_bin+=sigmaz_0_stat;
  sigmaz_0_et_bin+=sigmaz_0_et;


  psistarpsi_et_lby2_bin+=psistarpsi_et_lby2;
  psistarpsi_et_lby3_bin+=psistarpsi_et_lby3;
  psistarpsi_et_lby6_bin+=psistarpsi_et_lby6;
  psistarpsi_stat_lby2_bin+=psistarpsi_stat_lby2;
  psistarpsi_stat_lby3_bin+=psistarpsi_stat_lby3;
  psistarpsi_stat_lby6_bin+=psistarpsi_stat_lby6;
  tautau_et_lby2_bin+=tautau_et_lby2;
  tautau_et_lby3_bin+= tautau_et_lby3;
  tautau_et_lby6_bin+=tautau_et_lby6;
  tautau_stat_lby2_bin+=tautau_stat_lby2;
  tautau_stat_lby3_bin+=tautau_stat_lby3;
  tautau_stat_lby6_bin+=tautau_stat_lby6;
  psistarpsi_est_lby2_bin+=psistarpsi_est_lby2;
  psistarpsi_est_lby3_bin+=psistarpsi_est_lby3;
  psistarpsi_est_lby6_bin+= psistarpsi_est_lby6;
  psi2starpsi2_est_lby2_bin+=psi2starpsi2_est_lby2;
  psi2starpsi2_est_lby3_bin+=psi2starpsi2_est_lby3;
  psi2starpsi2_est_lby6_bin+=psi2starpsi2_est_lby6;
  psi3starpsi3_est_lby2_bin+=psi3starpsi3_est_lby2;
  psi3starpsi3_est_lby3_bin+=psi3starpsi3_est_lby3;
  psi3starpsi3_est_lby6_bin+=psi3starpsi3_est_lby6;
  tautau_est_lby2_bin+=tautau_est_lby2;
  tautau_est_lby3_bin+=tautau_est_lby3;
  tautau_est_lby6_bin+=tautau_est_lby6;
  psi3tau_est_lby2_bin+=psi3tau_est_lby2;
  psi3tau_est_lby3_bin+=psi3tau_est_lby3;
  psi3tau_est_lby6_bin+=psi3tau_est_lby6;
  phase_psi_bin+=phase_psi;
  repsi_bin+=repsi;
  impsi_bin=impsi;



  //incrementing av variables

  energy_av+=energy;
  sigmax_q_stat_av+=sigmax_q_stat;
  sigmaz_q_stat_av+= sigmaz_q_stat;
  sigmaz_q_et_av+=sigmaz_q_et;
  sigmaz_0_stat_av+=sigmaz_0_stat;
  sigmaz_0_et_av+=sigmaz_0_et;


  psistarpsi_et_lby2_av+=psistarpsi_et_lby2;
  psistarpsi_et_lby3_av+=psistarpsi_et_lby3;
  psistarpsi_et_lby6_av+=psistarpsi_et_lby6;
  psistarpsi_stat_lby2_av+=psistarpsi_stat_lby2;
  psistarpsi_stat_lby3_av+=psistarpsi_stat_lby3;
  psistarpsi_stat_lby6_av+=psistarpsi_stat_lby6;
  tautau_et_lby2_av+=tautau_et_lby2;
  tautau_et_lby3_av+= tautau_et_lby3;
  tautau_et_lby6_av+=tautau_et_lby6;
  tautau_stat_lby2_av+=tautau_stat_lby2;
  tautau_stat_lby3_av+=tautau_stat_lby3;
  tautau_stat_lby6_av+=tautau_stat_lby6;
  psistarpsi_est_lby2_av+=psistarpsi_est_lby2;
  psistarpsi_est_lby3_av+=psistarpsi_est_lby3;
  psistarpsi_est_lby6_av+= psistarpsi_est_lby6;
  psi2starpsi2_est_lby2_av+=psi2starpsi2_est_lby2;
  psi2starpsi2_est_lby3_av+=psi2starpsi2_est_lby3;
  psi2starpsi2_est_lby6_av+=psi2starpsi2_est_lby6;
  psi3starpsi3_est_lby2_av+=psi3starpsi3_est_lby2;
  psi3starpsi3_est_lby3_av+=psi3starpsi3_est_lby3;
  psi3starpsi3_est_lby6_av+=psi3starpsi3_est_lby6;
  tautau_est_lby2_av+=tautau_est_lby2;
  tautau_est_lby3_av+=tautau_est_lby3;
  tautau_est_lby6_av+=tautau_est_lby6;
  psi3tau_est_lby2_av+=psi3tau_est_lby2;
  psi3tau_est_lby3_av+=psi3tau_est_lby3;
  psi3tau_est_lby6_av+=psi3tau_est_lby6;
  phase_psi_av+=phase_psi;
  repsi_av+=repsi;
  impsi_av+=impsi;


  if(n_measure%5000==0)
  {
    cumf=fopen(cum_op_fname,"a");
    fprintf(cumf,"%li %.16le %.16le %.16le %.16le %.16le %.16le " ,
	n_measure,
	energy_av/((double)n_measure),
	sigmax_q_stat_av/((double)n_measure),
	sigmaz_q_stat_av/((double)n_measure),
	sigmaz_q_et_av/((double)n_measure),
	sigmaz_0_stat_av/((double)n_measure),
	sigmaz_0_et_av/((double)n_measure));
    fprintf(cumf,"\n");
    fclose(cumf);

    cumf=fopen(cum_op_local_fname,"a");
    fprintf(cumf,"%li  %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le fake %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le ",
	n_measure,
	psistarpsi_et_lby2_av/((double)n_measure),
	psistarpsi_et_lby3_av/((double)n_measure),
	psistarpsi_et_lby6_av/((double)n_measure),
	psistarpsi_stat_lby2_av/((double)n_measure),
	psistarpsi_stat_lby3_av/((double)n_measure),
	psistarpsi_stat_lby6_av/((double)n_measure),
	tautau_et_lby2_av/((double)n_measure),
	tautau_et_lby3_av/((double)n_measure),
	tautau_et_lby6_av/((double)n_measure),
	tautau_stat_lby2_av/((double)n_measure),
	tautau_stat_lby3_av/((double)n_measure),
	tautau_stat_lby6_av/((double)n_measure),
	psistarpsi_est_lby2_av/((double)n_measure),
	psistarpsi_est_lby3_av/((double)n_measure),
	psistarpsi_est_lby6_av/((double)n_measure),
	psi2starpsi2_est_lby2_av/((double)n_measure),
	psi2starpsi2_est_lby3_av/((double)n_measure),
	psi2starpsi2_est_lby6_av/((double)n_measure),
	psi3starpsi3_est_lby2_av/((double)n_measure),
	psi3starpsi3_est_lby3_av/((double)n_measure),
	psi3starpsi3_est_lby6_av/((double)n_measure),
	tautau_est_lby2_av/((double)n_measure),
	tautau_est_lby3_av/((double)n_measure),
	tautau_est_lby6_av/((double)n_measure),
	psi3tau_est_lby2_av/((double)n_measure),
	psi3tau_est_lby3_av/((double)n_measure),
	psi3tau_est_lby6_av/((double)n_measure),
	phase_psi_av/((double)n_measure),
	repsi_av/((double)n_measure),
	impsi_av/((double)n_measure));
    fprintf(cumf,"\n");
    fclose(cumf);



    histf=fopen(hist_fname,"w");
    for(i=0;i<63;i++) {
      fprintf(histf,"%i %le\n",i,(double)hist_phase[i]/(double)n_measure);
    }


  }//n_measure%5000==0

  if(n_measure%(max_mc_step/6)==0){


    clustf=fopen(clusterst_fname,"w"); 
    hist_clust_size[N_CLUST_HIST_BINS-1]+=hist_clust_size[N_CLUST_HIST_BINS]; 
    for(i=0;i<(N_CLUST_HIST_BINS);i++) 
      if(hist_clust_size[i]!=0){ 
	fprintf(clustf,"%d\t%.15e\n",i,(((double)hist_clust_size[i])/((double)num_cluster))); 
      }
    fclose(clustf); 
  }

if(num_call_bin % binsize == 0)
{
  completed_bins++;
  binf=fopen(bin_op_fname,"a");
  fprintf(binf,"%li %li %li %.16le %.16le %.16le %.16le %.16le %.16le " ,
      completed_bins,
      n_measure,
      num_call_bin,
      energy_bin/((double)binsize),
      sigmax_q_stat_bin/((double)binsize),
      sigmaz_q_stat_bin/((double)binsize),
      sigmaz_q_et_bin/((double)binsize),
      sigmaz_0_stat_bin/((double)binsize),
      sigmaz_0_et_bin/((double)binsize));
  fprintf(binf,"\n");
  fclose(binf);

  binf=fopen(bin_op_local_fname,"a");
  fprintf(binf,"%li %li %li %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le ",
      completed_bins,
      n_measure,
      num_call_bin,
      psistarpsi_et_lby2_bin/((double)binsize),
      psistarpsi_et_lby3_bin/((double)binsize),
      psistarpsi_et_lby6_bin/((double)binsize),
      psistarpsi_stat_lby2_bin/((double)binsize),
      psistarpsi_stat_lby3_bin/((double)binsize),
      psistarpsi_stat_lby6_bin/((double)binsize),
      tautau_et_lby2_bin/((double)binsize),
      tautau_et_lby3_bin/((double)binsize),
      tautau_et_lby6_bin/((double)binsize),
      tautau_stat_lby2_bin/((double)binsize),
      tautau_stat_lby3_bin/((double)binsize),
      tautau_stat_lby6_bin/((double)binsize),
      psistarpsi_est_lby2_bin/((double)binsize),
      psistarpsi_est_lby3_bin/((double)binsize),
      psistarpsi_est_lby6_bin/((double)binsize),
      psi2starpsi2_est_lby2_bin/((double)binsize),
      psi2starpsi2_est_lby3_bin/((double)binsize),
      psi2starpsi2_est_lby6_bin/((double)binsize),
      psi3starpsi3_est_lby2_bin/((double)binsize),
      psi3starpsi3_est_lby3_bin/((double)binsize),
      psi3starpsi3_est_lby6_bin/((double)binsize),
      tautau_est_lby2_bin/((double)binsize),
      tautau_est_lby3_bin/((double)binsize),
      tautau_est_lby6_bin/((double)binsize),
      psi3tau_est_lby2_bin/((double)binsize),
      psi3tau_est_lby3_bin/((double)binsize),
      psi3tau_est_lby6_bin/((double)binsize),
      phase_psi_bin/((double)binsize),
      repsi_bin/((double)binsize),
      impsi_bin/((double)binsize));
  fprintf(binf,"\n");
  fclose(binf);


}

}
