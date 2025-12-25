#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"global.h"
#include<unistd.h>
#include"time.h"
#include"mt19937ar.h"
void makelinks(int);
void freelinks();
void do_diag(int);
void do_clust(void);
void set_const_val(void);
void make_lattice(void);
void do_initialisation(void);
void do_metropolis_ising(void);
void measure(void);
void fat_update(int);
void free_fat_graph();
void fat_update(int typ);
void make_fat_graph();
void canonical();

int main()
{
    tmake=tclust=0;
    tstitch=0;
    unsigned long mc_step;
    int i,k,counter,xsector,ysector;
    long int hist_xsector[LX+1],hist_ysector[LY+1];	
    set_const_val();
    make_lattice();
    do_initialisation();
    mc_step=0;
    int ifwup;
    printf("tot=%lu\n",totalmcsteps);
    time_t t0,t1;

    
    for(mc_step=0; mc_step<(totalmcsteps);mc_step++){
        ifwup=(mc_step<max_w_step);
        if(mc_step%5000==0)
            printf("%lu \n",mc_step);



        for(k=0;k<3;k++){
            do_diag(ifwup);
            makelinks(k);
            wolffsteps=(n_niop-n_triagop);
            canonical();
            freelinks();

            //makelinks(k);
            make_fat_graph();

            fat_update(k);
            free_fat_graph();
            //do_clust();
            //freelinks();
        }
            if(!ifwup)
                measure();

    }	
    printf("done steps, make=%f, clust=%f, stitch=%f\n",tmake,tclust,tstitch);

    printf("done %li steps\n",totalmcsteps);
    free(opstr);
    free(divider);
    free(minority);
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



