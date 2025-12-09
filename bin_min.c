#include<stdio.h>
#include "global.h"
void bin(double energy,double sigmax_q_stat,double sigmaz_q_stat,double sigmaz_q_et,double sigmaz_0_stat,double sigmaz_0_et){



  FILE *binf, *cumf;

  static long num_call_bin=0;

  static double  energy_bin;
  static double sigmax_q_stat_bin;
  static double sigmaz_q_stat_bin;
  static double sigmaz_q_et_bin;
  static double sigmaz_0_stat_bin;
  static double sigmaz_0_et_bin;



  //initialise bin variables

  if(!(num_call_bin % binsize))
  {
    energy_bin=0.0;
   sigmax_q_stat_bin=0.0;
   sigmaz_q_stat_bin=0.0;
   sigmaz_q_et_bin=0.0;
   sigmaz_0_stat_bin=0.0;
   sigmaz_0_et_bin=0.0;



  }

  num_call_bin++;	

  //incrementing bin variables

   energy_bin+=energy;
  sigmax_q_stat_bin+=sigmax_q_stat;
   sigmaz_q_stat_bin+= sigmaz_q_stat;
   sigmaz_q_et_bin+=sigmaz_q_et;
   sigmaz_0_stat_bin+=sigmaz_0_stat;
   sigmaz_0_et_bin+=sigmaz_0_et;




  //incrementing av variables
  
   energy_av+=energy;
   sigmax_q_stat_av+=sigmax_q_stat;
   sigmaz_q_stat_av+= sigmaz_q_stat;
   sigmaz_q_et_av+=sigmaz_q_et;
   sigmaz_0_stat_av+=sigmaz_0_stat;
   sigmaz_0_et_av+=sigmaz_0_et;



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
    



  }//n_measure%5000==0


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


  }

}
