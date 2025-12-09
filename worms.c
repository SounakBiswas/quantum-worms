#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
int get_sdl_from_sdv(int sdv, int sdvnbr);
int get_sl_at_sdl(int sdl);
int get_sdl_at_sl(int sl);
int get_sl_from_sv(int sv,int svnbr);
//graph data structures
int *sv_at_v;
int *v_at_sv;
int *dv_at_sdv;
int *sdv_at_dv;
int *l_at_sl;
int *dl_at_sdl;
int *l_at_dl;
int *p_at_sp;
int *l2v[2];
int *dl2dv[2];
int *dv2v[3];
int *dv_set;
int *dv_nbrctr;
int **dv_nbr; //nbr table for dual vertices
int **dv2dl; //site to link for dual vertices //int **dv2dl;
int *fsp;
int *ncr_at_l;
int *dimer;
void create_graph(){
    // dual graphs and links have prefix d
    // spatial and spatial-dual are s and sd
    int site;
    int vctr=0;
    int lctr=0;
    int dlctr=0;
    int dvctr=0;
    int pctr=0;
    int *svert;
    int op;
    int s1,s2,s3;
    int dsite;
    int sdv;
    int *dv_typ; //wt type
    //add direct graph
    for (site=0; site<nsites; site++){
        sv_at_v[vctr]=site;
        v_at_sv[site]=vctr;
        fsp[vctr]=sigma[site];
        vctr++;
    }
    int sv;
    for (sv=0; sv<nsites; sv++){
        for(int i=0; i<3; i++){
            int svnbr = nbr[site][i];
            int vnbr = v_at_sv[svnbr];
            int sl = get_sl_from_sv(sv,svnbr);
            l_at_sl[sl] = lctr;
            l2v[0][lctr] = v_at_sv[sv];
            l2v[1][lctr] = v_at_sv[svnbr];
            ncr_at_l[lctr]=0;
            lctr++;
        }
    }
    // add dual graph
    for (sdv=0; sdv<ndsites; sdv++){
        sdv_at_dv[dvctr]=sdv;
        dv_at_sdv[dsite]=dvctr;
        dv_typ[dvctr]=0;
        dvctr++;
    }
    for (sdv=0; dsite<ndsites; sdv++){
        for(int psite=0; psite<npsites ;psite++){
            sv = triangle[dsite][psite];
            dv2v[psite][dvctr]= v_at_sv[sv];
        }
        for (int i=0; i<nsdnbrs; i++){
            int sdvnbr=dnbr[sdv][i];
            int sdl = get_sdl_from_sdv(sdv,sdvnbr);
            int sl= get_sl_at_sdl(sdl);
            int l=l_at_sl[sl];
            if(ncr_at_l[l]==0){
                int dv=dv_at_sdv[sdv];
                int dvnbr=dv_at_sdv[sdvnbr];
                dl2dv[0][dlctr] = dv;
                dl2dv[1][dlctr] = dvnbr;
                dv2dl[dv][dv_nbrctr[dv]++]=dlctr;
                dv2dl[dvnbr][dv_nbrctr[dvnbr]++]=dlctr;
                l_at_dl[dlctr]=l;
                dlctr++;
                ncr_at_l[l]=1;
            }
            else ncr_at_l[l]=-1;
        }
    }

    for (int op_pos=0; op_pos <opstr_l; op_pos++) {
        op = opstr[op_pos];
        //new segment
        if(op < ntriangles){
            int splaq=op;
            int sdv=splaq;
            int dv = dv_at_sdv[sdv];
            dv_typ[dv]=1;


            if(dv_set[sdv]==0){

                //set direct graph

                //set dual graph (**crucial**)
                sdv_at_dv[dvctr]=sdv;
                dv_at_sdv[sdv]=dvctr;
                dvctr++;
                for(int tsite=0; tsite<3; tsite++){
                    site = triangle[op][tsite];
                    dv2v[tsite][dvctr]= v_at_sv[site];
                }
                //add 3 links 
                for (int i=0; i<nsdnbrs; i++){
                    int sdvnbr=dnbr[sdv][i];
                    int sdl = get_sdl_from_sdv(sdv,sdvnbr);
                    int sl= get_sl_at_sdl(sdl);
                    int l=l_at_sl[sl];
                    if(ncr_at_l[l]==0){
                        int dv=dv_at_sdv[sdv];
                        int dvnbr=dv_at_sdv[sdvnbr];
                        dl2dv[0][dlctr] = dv;
                        dl2dv[1][dlctr] = dvnbr;
                        dv2dl[dv][dvdlctr[dv]++]=dlctr;
                        dv2dl[dvnbr][dvdlctr[dvnbr]++]=dlctr;
                        l_at_dl[dlctr]=l;
                        dlctr++;
                        ncr_at_l[l]=1;
                    }
                    else ncr_at_l[l]=-1;
                }

            }


        }
        if (op >= ntriangles){
            sv= (op-ntriangles)%nsites;
            if(op>ndiagops){
                sigma[sv]*=-1;
            }

            // add new vertex
            //modify_direct_graph(){
            //}
            v_at_sv[sv]=vctr;
            fsp[vctr] = sigma[sv];
            vctr++;
            //add links for new vertex (crucial to check for hyperedges
            for(int i=0; i<6; i++){
                int svnbr = nbr[sv][i];
                int vnbr = v_at_sv[svnbr];
                int sl=get_sl_from_sv(sv,svnbr);
                l_at_sl[sl] = lctr;
                if(i<3){
                    l2v[0][lctr] = v_at_sv[sv];
                    l2v[1][lctr] = v_at_sv[svnbr];
                }
                else {
                    l2v[0][lctr] = v_at_sv[sv];
                    l2v[1][lctr] = v_at_sv[svnbr];
                }
                ncr_at_l[lctr]=0;
                lctr++;
            }
            

            //dual graph vertices
            //modify_dual_graph(sv);
            //add 6 vertices
            for (int i=0; i<nplaqspersite; i++){
              sdv=sdv_at_sv[sv][i];
              sdv_at_dv[dvctr]=sdv;
              dvctr++;
              for(int tsite=0; tsite<3; tsite++){
                  sv = triangle[sdv][tsite];
                  dv2v[tsite][dvctr]= v_at_sv[site];
              }
              //***add triangles if using, dv2v

            }
            for (int i=0; i<nplaqspersite; i++){
              sdv=sdv_at_sv[sv][i];
              //uptriangles
              if(sdv%2==0){
                int sdv0=dnbr[sdv][0];
                //add_link(sdv,sdv0)
                int sdv1=dnbr[sdv][1];
                //add_link(sdv,sdv1)
                int sdv2=dnbr[sdv][2];
                //add_link(sdv,sdv2)

              }
              else{
                int sdv2=dnbr[sdv][2];
                //add_link(sdv2,sdv)
              }
             


            }




        }
    }

}
void add_link( int sdv0, int sdv1,int *dlctr){
    int sdl = get_sdl_from_sdv(sdv0, sdv1);
    int sl= get_sl_at_sdl(sdl);
    int l=l_at_sl[sl];
    //check if it corresponds to a hyperedge
    if(ncr_at_l[l]==0){
        int dv0=dv_at_sdv[sdv0];
        int dv1=dv_at_sdv[sdv1];
        dl2dv[0][*dlctr] = dv0;
        dl2dv[1][*dlctr] = dv1;
        dv2dl[dv0][dv_nbrctr[dv0]]=(*dlctr);
        dv2dl[dv1][dv_nbrctr[dv1]]=(*dlctr);
        dv_nbr[dv0][dv_nbrctr[dv0]]=dv1;
        dv_nbr[dv1][dv_nbrctr[dv1]]=dv0;
        dv_nbrctr[dv0]++;
        dv_nbrctr[dv1]++;
        l_at_dl[*dlctr]=l;
        (*dlctr)++;
        ncr_at_l[l]=1;
    }
    else ncr_at_l[l]=-1;
}


