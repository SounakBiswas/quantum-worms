#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"worms.h"
//graph data structures
int *sv_at_v;
int *v_at_sv;
int *dv_at_sdv;
int *sdv_at_dv;
int *l_at_sl;
int *dl_at_sdl;
int *l_at_dl;
int *dl_at_l;
int *p_at_sp;
int *v_at_l[2];
int *dv_at_dl[2];
int *dv2v[3];
int *dv_set;
int *dv_nbrctr;
int *v_nbrctr;
int **l_at_v;
int *dl_at_dv[3]; //site to link for dual vertices //int **dl_at_dv;
int *fsp;
int *if_hyperedge;
int *dimer;
int *dv_typ; //wt type
int vctr;
int lctr;
int dlctr;
int dvctr;
int *wx_mark;
int *wy_mark;
void add_l_to_v(int v,int lctr){
  v_nbrctr[v]++;
  l_at_v[v]=realloc(l_at_v[v],sizeof(int)*v_nbrctr[v]);
}
void create_graph(){
    // dual graphs and links have prefix d
    // spatial and spatial-dual are s and sd
    int site;
    int vctr=0;
    lctr=0;
    dlctr=0;
    dvctr=0;
    int pctr=0;
    int op;
    int s1,s2,s3;
    int dsite;
    int sdv;
    int v0,v1;
    //add direct graph
    for (site=0; site<nsites; site++){
        sv_at_v[vctr]=site;
        v_at_sv[site]=vctr;
        fsp[vctr]=sigma[site];
        dv_typ[vctr]=0;
        vctr++;
    }
    int sv;
    for (sv=0; sv<nsites; sv++){
        for(int i=0; i<3; i++){
            int svnbr = nbr[site][i];
            int vnbr = v_at_sv[svnbr];
             
            int sl = get_sl_from_sv(sv,svnbr);
            l_at_sl[sl] = lctr;
            wx_mark[lctr]=wx_smark[sl];
            wy_mark[lctr]=wy_smark[sl];
            v0=v_at_sv[sv];
            v1=v_at_sv[svnbr];
            dimer[lctr]= (fsp[v0]==fsp[v1]) ? 1 : 0;
            if_hyperedge[lctr]=0;
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
        
        for (int i=0; i<nplaqspersite; i++){
          sdv=sdv_at_sv[sv][i];
          //add dual links for uptriangles
          if(sdv%2==0){
              int sdv0=dnbr[sdv][0];
              add_dual_link(sdv,sdv0,&dlctr);
              int sdv1=dnbr[sdv][1];
              add_dual_link(sdv,sdv1,&dlctr);
              int sdv2=dnbr[sdv][2];
              add_dual_link(sdv,sdv2,&dlctr);
          }

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




        }
        if (op >= ntriangles){
            sv= (op-ntriangles)%nsites;
            if(op>ndiagops){
                sigma[sv]*=-1;
            }


            v_at_sv[sv]=vctr;
            fsp[vctr] = sigma[sv];
            v0=v_at_sv[sv];
            dv_typ[vctr]=0;
            vctr++;
            //add links for new vertex (crucial to check for hyperedges
            for(int i=0; i<6; i++){
                int svnbr = nbr[sv][i];
                int vnbr = v_at_sv[svnbr];
                int sl=get_sl_from_sv(sv,svnbr);
                v1=v_at_sv[svnbr];
                l_at_sl[sl] = lctr;
                if(i<3){
                    v_at_l[0][lctr] = v0;
                    v_at_l[1][lctr] = v1;
                }
                else {
                    v_at_l[0][lctr] = v0;
                    v_at_l[1][lctr] = v1;
                }
                l_at_v[v0][v_nbrctr[v0]++]=lctr;
                l_at_v[v1][v_nbrctr[v1]++]=lctr;
                dimer[lctr]= (fsp[v0]==fsp[v1]) ? 1 : 0;
                if_hyperedge[lctr]=0;
                lctr++;
            }
            //dual graph vertices
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
                  add_dual_link(sdv,sdv0,&dlctr);
                  int sdv1=dnbr[sdv][1];
                  add_dual_link(sdv,sdv1,&dlctr);
                  int sdv2=dnbr[sdv][2];
                  add_dual_link(sdv,sdv2,&dlctr);
              }
              else{
                  int sdv2=dnbr[sdv][2];
                  add_dual_link(sdv2,sdv,&dlctr);
              }

            }
        }
    }

}
void add_dual_link( int sdv0, int sdv1,int *dlctr){
    int sdl = get_sdl_from_sdv(sdv0, sdv1);
    int sl= get_sl_at_sdl(sdl);
    int l=l_at_sl[sl];
    //check if it corresponds to a hyperedge
    if(if_hyperedge[l]==-1){
        int dv0=dv_at_sdv[sdv0];
        int dv1=dv_at_sdv[sdv1];
        dv_at_dl[0][*dlctr] = dv0;
        dv_at_dl[1][*dlctr] = dv1;
        dl_at_dv[dv0][dv_nbrctr[dv0]]=(*dlctr);
        dl_at_dv[dv1][dv_nbrctr[dv1]]=(*dlctr);
        dv_nbr[dv0][dv_nbrctr[dv0]]=dv1;
        dv_nbr[dv1][dv_nbrctr[dv1]]=dv0;
        dv_nbrctr[dv0]++;
        dv_nbrctr[dv1]++;
        l_at_dl[*dlctr]=l;
        dl_at_l[l]= (*dlctr);
        (*dlctr)++;
    }
    if_hyperedge[l]++;
}


