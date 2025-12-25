#include<stdio.h>
#include<math.h>
#include"global.h"

void set_const_val(void)
{
  char file_path[120];
  int ext;

 root3by2=sqrt(3)*0.5;
  pi=M_PI;
  piby2=pi/2;
 
  //parameters
 
  //*************************
  max_mc_step = MAX_MC_STEP;
  max_w_step = MAX_W_STEP;
  totalmcsteps = max_mc_step+max_w_step;
  printf("%lu %lu %lu \n", max_mc_step, max_w_step,totalmcsteps);
  binsize = BINSIZE;

  tau_max = TAU_MAX;
  //***************************

  nsites=NSITES;
  nbonds=NBONDS;
  ntrianglesandsites=NTRIANGLESANDSITES;
  ntriangles=NTRIANGLES;
  ndvperv=NDVPERV;
  nplaqspersite=6;
  ndiagops=ntrianglesandsites;
  //dual lattice
  ndsites=NDSITES;
  nsdl=NDLINKS;
  nsdnbrs=3;
  lx=LX;ly=LY;
  initialstate=INITIALSTATE;
  temp=TEMP;
  beta=1.0/temp;
  h=H;
  jising=JISING; 

  


  //filenames
  ext=45;
  sprintf(file_path,"./outfiles/");

  sprintf(bin_op_fname,"%sbinopL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(cum_op_fname,"%scumopL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(bin_op_local_fname,"%sbinlocalL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(cum_op_local_fname,"%scumlocalL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(auto_fname,"%sautocorrL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(hist_fname,"%shistL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(clusterst_fname,"%sclusterstL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  sprintf(save_fname,"%ssavL%ij%.4lft%.4lfh%.4lf",file_path,lx,jising,temp,H);

  Sin[0]=0;Sin[1]=root3by2;Sin[2]=-root3by2;
  Cos[0]=1;Cos[1]=-0.5;Cos[2]=-0.5;

  return;
}

