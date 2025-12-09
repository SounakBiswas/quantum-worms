#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"global.h"
#include<unistd.h>
#include"time.h"
#include"mt19937ar.h"
void makelinks(int);
void do_diag(void);
void do_diag_withadjust(void);
void do_clust(void);
void set_const_val(void);
void make_lattice(void);
void do_initialisation(void);
void do_metropolis_ising(void);
void measure(void);

int main()
{
  double t,t1,t2; t=0; // for time measurements
  unsigned long mc_step;
  char stoppingname[256];
  FILE *stop;
  int i,k,counter,xsector,ysector;
  long int hist_xsector[LX+1],hist_ysector[LY+1];	
  set_const_val();
  make_lattice();
  sprintf(stoppingname,"stopfileL%ij%.4lft%.4lfh%.4lf",lx,jising,temp,H);
  /* Function below from unistd.h checks if file with name saved in save_fname exists. If it exists, then it is read to get all initialization information. Else it is assumed that the run is begining and all variables are initialized. So delete sav files from directory outfiles before running unless the run is a continuation of a previous monte-carlo chain. */

    do_initialisation();
    mc_step=0;

  if(mc_step > 0) 
    printf("run restarted after %lu steps\n",mc_step);
  //printf("%d\t%d\n",n_niop,n_triagop);
  for(;mc_step<(totalmcsteps);mc_step++)
  {

 if(mc_step < max_w_step){
   for(k=0;k<3;k++){
    do_diag_withadjust();
    makelinks(k);
    t1=clock();
    do_clust();
    t2=clock();
    t+=(double)(t2-t1)/CLOCKS_PER_SEC;
   }
 } else {
   for(k=0;k<3;k++){
    do_diag();
    makelinks(k);
    t1=clock();
    do_clust();
    t2=clock();
    t+=(double)(t2-t1)/CLOCKS_PER_SEC;
   }
 }
   //printf("%lu\n",mc_step);

    if(mc_step==(max_w_step-1)){ 
      printf("warmup completed with %lu steps\n",mc_step+1);
      hist_clust_size=(int*)malloc((N_CLUST_HIST_BINS+1)*sizeof(int));
      for(i=0;i<=N_CLUST_HIST_BINS;i++)
	hist_clust_size[i]=0;
    take_stat_bool=1;
    }
    if(mc_step > (max_w_step-1))
    {
      //	  checker();
     if((mc_step - max_w_step+1)%2 == 0)
       measure();

    }


  }	

 printf("done %li steps\n",totalmcsteps);

  if(totalmcsteps==max_w_step+max_mc_step)
  {
    stop=fopen(stoppingname,"w");
    fprintf(stop,"stopped ran=%lf",genrand_real2());
    fclose(stop);


  }
  free(opstr);
  /*write cluster statisitcs*/
  /*************************************/
  free(hist_clust_size) ; 
  /* write time data for time */
//  FILE* ti;
 // ti=fopen("speed","a");
 // fprintf(ti,"%d\t%f\n",nsites,t);
 //fclose(ti);
  /* ******************* */
  return 0;

}



