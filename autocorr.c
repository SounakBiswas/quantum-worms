#include<stdio.h>
#include"global.h"
#include"math.h"

#define min(a,b) ((a<b)?a:b)

void auto_corr(double op1,double op2,double op3,double op4,double op5)
{

   FILE *autof;

   static long int n_autocorr=0;
   static int counter[TAU_MAX]={0};

   static double stored_op1[TAU_MAX];
   static double stored_op2[TAU_MAX];
   static double stored_op3[TAU_MAX];
   static double stored_op4[TAU_MAX];
   static double stored_op5[TAU_MAX];

   static double ac_op1[TAU_MAX]={0};
   static double ac_op2[TAU_MAX]={0};
   static double ac_op3[TAU_MAX]={0};
   static double ac_op4[TAU_MAX]={0};
   static double ac_op5[TAU_MAX]={0};

   static double acerr_op1[TAU_MAX]={0};
   static double acerr_op2[TAU_MAX]={0};
   static double acerr_op3[TAU_MAX]={0};
   static double acerr_op4[TAU_MAX]={0};
   static double acerr_op5[TAU_MAX]={0};

   static double av_op1, av_op2, av_op3, av_op4, av_op5;

   //local variables
   int t,t1,tau;

   t=n_autocorr%tau_max;

   stored_op1[t]=op1;
   stored_op2[t]=op2;
   stored_op3[t]=op3;
   stored_op4[t]=op4;
   stored_op5[t]=op5;

   av_op1 += op1;
   av_op2 += op2;
   av_op3 += op3;
   av_op4 += op4;
   av_op5 += op5;


   for(tau=0; tau<min(n_autocorr,tau_max); tau++) {
      t1=(n_autocorr-tau)%tau_max;
      counter[tau]++;

      ac_op1[tau] += stored_op1[t]*stored_op1[t1];
      ac_op2[tau] += stored_op2[t]*stored_op2[t1];
      ac_op3[tau] += stored_op3[t]*stored_op3[t1];
      ac_op4[tau] += stored_op4[t]*stored_op4[t1];
      ac_op5[tau] += stored_op5[t]*stored_op5[t1];

      acerr_op1[tau] += (stored_op1[t]*stored_op1[t1])*(stored_op1[t]*stored_op1[t1]);
      acerr_op2[tau] += (stored_op2[t]*stored_op2[t1])*(stored_op2[t]*stored_op2[t1]);
      acerr_op3[tau] += (stored_op3[t]*stored_op3[t1])*(stored_op3[t]*stored_op3[t1]);
      acerr_op4[tau] += (stored_op4[t]*stored_op4[t1])*(stored_op4[t]*stored_op4[t1]);
      acerr_op5[tau] += (stored_op5[t]*stored_op5[t1])*(stored_op5[t]*stored_op5[t1]);
   }
   n_autocorr++;	

   if(n_autocorr == (max_mc_step/2) ) {
      for(tau=0;tau<tau_max;tau++) {
         ac_op1[tau] /= (double)counter[tau];
         ac_op2[tau] /= (double)counter[tau];
         ac_op3[tau] /= (double)counter[tau];
         ac_op4[tau] /= (double)counter[tau];
         ac_op5[tau] /= (double)counter[tau];

         acerr_op1[tau] /= (double)counter[tau];
         acerr_op2[tau] /= (double)counter[tau];
         acerr_op3[tau] /= (double)counter[tau];
         acerr_op4[tau] /= (double)counter[tau];
         acerr_op5[tau] /= (double)counter[tau];

         acerr_op1[tau] -= ac_op1[tau]*ac_op1[tau];
         acerr_op2[tau] -= ac_op2[tau]*ac_op2[tau];
         acerr_op3[tau] -= ac_op3[tau]*ac_op3[tau];
         acerr_op4[tau] -= ac_op4[tau]*ac_op4[tau];
         acerr_op5[tau] -= ac_op5[tau]*ac_op5[tau];

         acerr_op1[tau] = sqrt(acerr_op1[tau])/sqrt((double)counter[tau]-1.0);
         acerr_op2[tau] = sqrt(acerr_op2[tau])/sqrt((double)counter[tau]-1.0);
         acerr_op3[tau] = sqrt(acerr_op3[tau])/sqrt((double)counter[tau]-1.0);
         acerr_op4[tau] = sqrt(acerr_op4[tau])/sqrt((double)counter[tau]-1.0);
         acerr_op5[tau] = sqrt(acerr_op5[tau])/sqrt((double)counter[tau]-1.0);

         ac_op1[tau] -= av_op1*av_op1/((double)(n_autocorr*n_autocorr));
         ac_op2[tau] -= av_op2*av_op2/((double)(n_autocorr*n_autocorr));
         ac_op3[tau] -= av_op3*av_op3/((double)(n_autocorr*n_autocorr));
         ac_op4[tau] -= av_op4*av_op4/((double)(n_autocorr*n_autocorr));
         ac_op5[tau] -= av_op5*av_op5/((double)(n_autocorr*n_autocorr));
      }

      autof=fopen(auto_fname,"w");

      for(tau=0;tau<tau_max;tau++) {
         fprintf(autof,"%d %le %le %le %le %le %le %le %le %le %le\n",tau,ac_op1[tau]/ac_op1[0],acerr_op1[tau]/ac_op1[0],ac_op2[tau]/ac_op2[0],acerr_op2[tau]/ac_op2[0],ac_op3[tau]/ac_op3[0],acerr_op3[tau]/ac_op3[0],ac_op4[tau]/ac_op4[0],acerr_op4[tau]/ac_op4[0],ac_op5[tau]/ac_op5[0],acerr_op5[tau]/ac_op5[0]); }
      fclose(autof);
   }
}
