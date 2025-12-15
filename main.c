#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"global.h"
#include<unistd.h>
#include"time.h"
#include"mt19937ar.h"
void makelinks(int);
void do_diag(int);
void do_clust(void);
void set_const_val(void);
void make_lattice(void);
void do_initialisation(void);
void do_metropolis_ising(void);
void measure(void);
void worm_update();

int main()
{
    double t,t1,t2; t=0; // for time measurements
    unsigned long mc_step;
    int i,k,counter,xsector,ysector;
    long int hist_xsector[LX+1],hist_ysector[LY+1];	
    set_const_val();
    make_lattice();
    do_initialisation();
    mc_step=0;
    int ifwup;

    for(mc_step=0; mc_step<(totalmcsteps);mc_step++){
        ifwup=(mc_step<max_w_step);



        for(k=0;k<3;k++){
            do_diag(ifwup);
            makelinks(k);
            t1=clock();
            do_clust();
            t2=clock();
            t+=(double)(t2-t1)/CLOCKS_PER_SEC;
        }

    }	
    printf("done steps\n");
    worm_update();

    printf("done %li steps\n",totalmcsteps);
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



