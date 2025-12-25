#ifndef GLOBAL_H
#define GLOBAL_H


//////parameters:
#ifdef DEBUG
    #define DEBUG_PRINT(...) printf("DEBUG: " __VA_ARGS__)
#else
    #define DEBUG_PRINT(...) do {} while (0)
#endif
#ifdef DEBUG
    #define DEBUG_GC() getchar()
#else
    #define DEBUG_GC() do {} while (0)
#endif


#define MAX_MC_STEP 200000
#define MAX_W_STEP 20000
#define BINSIZE 100
#define LX 9
#define LY 9
#define NSITES (LX*LY)
#define NTRIANGLES (2*LX*LY)
#define NBONDS (3*LX*LY)
#define NTRIANGLESANDSITES (3*LX*LY)
#define TEMP 1.0000
#define H 1.0000
#define NDVPERV 6
#define JISING 1.0000// Ising exchange
#define INITIALSTATE 0 //initial state 2: staggered, initialstate 1: columnar ; initialstate=0: all sigma random;
#define SEED 21
#define INIT_OPSTR_L 30
#define NDSITES 2*NSITES
#define NDLINKS 3*NSITES
#define TAU_MAX 1000
#define NLINKS NBONDS
#define N_CLUST_HIST_BINS 100
///////variables:
///dual graph
int sdv_at_sv[NSITES][6];
int sl_at_sdl[NDLINKS];
int wx_smark[NLINKS], wy_smark[NLINKS];
int dnbr[2*NSITES][3]; //gives a,b,c site numbers of given triangle
int ndvperv;
double tmake,tclust,tstitch;
                       //

//constants
double jising;
int nsites,nbonds,lx,ly,initialstate,ntriangles,ntrianglesandsites; 
int nsdnbrs, ndsites,npsites, nsdl;
int ndiagops;
int nplaqspersite;
double beta,h,temp;
double root3by2;
double piby2;
double pi;
int sublat[NSITES];
int triangle[2*NSITES][3]; //gives a,b,c site numbers of given triangle
int ptolink[2*NSITES][3];
                           //
//inputoutput
char auto_fname[256],bin_op_local_fname[256],cum_op_local_fname[256],bin_op_fname[256],cum_op_fname[256],save_fname[256],hist_fname[256],clusterst_fname[256];
//config and geometry
int sigma[NSITES];
int nbr[NSITES][6];
int site2[NSITES],site3[NSITES],pluslby2[NSITES],pluslby3[NSITES],pluslby6[NSITES],minuslby6[NSITES],minuslby2[NSITES],minuslby3[NSITES];
//Monte Carlo

//diagonal update 
//int opstr[opstr_l];
int opstr_l;
int new_opstr_l;
int *opstr;
int *divider;
int n_niop; //number of non identity operators in the seq
int n_triagop; // number of kedar triangles.
int n_flipop;
//clusterupdate
int legs; // nop of legs in the linked list    
int *op_pos;
int *leg_pos;
int *link_list;
int first[NSITES];
int last[NSITES];
int *hist_clust_size;
int take_stat_bool;
int num_cluster;


//measurement
unsigned long n_measure,completed_bins;
unsigned long totalmcsteps,binsize;
double Sin[3],Cos[3];
//*************************************************
unsigned long max_mc_step,max_w_step,numberofbins;
//*************************************************
unsigned long tau_max;
double energy_av;
double sigmaz_q_stat_av;
double sigmaz_q_et_av;
double sigmaz_0_et_av;
double sigmaz_0_stat_av;
double sigmax_q_stat_av;


double tautau_et_lby2_av,tautau_et_lby3_av,tautau_et_lby6_av;
double tautau_stat_lby2_av,tautau_stat_lby3_av,tautau_stat_lby6_av;
double tautau_est_lby2_av,tautau_est_lby3_av,tautau_est_lby6_av;
 

double psistarpsi_et_lby2_av,psistarpsi_et_lby3_av,psistarpsi_et_lby6_av;
double psistarpsi_stat_lby2_av,psistarpsi_stat_lby3_av,psistarpsi_stat_lby6_av;
double psistarpsi_est_lby2_av,psistarpsi_est_lby3_av,psistarpsi_est_lby6_av;

double psi2starpsi2_est_lby2_av,psi2starpsi2_est_lby3_av,psi2starpsi2_est_lby6_av;
double psi3starpsi3_est_lby2_av,psi3starpsi3_est_lby3_av,psi3starpsi3_est_lby6_av;
double psi3tau_est_lby2_av,psi3tau_est_lby3_av,psi3tau_est_lby6_av;

double phase_psi_av;
double impsi_av,repsi_av;



double hist_phase[63];


#endif //GLOBAL_H
