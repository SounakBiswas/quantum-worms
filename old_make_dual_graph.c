#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"worms.h"
//graph data structures
int *first_v;
int *first_dv;
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
int **dl_at_dv; //site to link for dual vertices //int **dl_at_dv;
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
int get_partner(int v1,int link);

void add_l_to_v(int v,int lctr){
    v_nbrctr[v]++;
    l_at_v[v]=realloc(l_at_v[v],sizeof(int)*v_nbrctr[v]);
}
void add_dl_to_dv(int dv,int dlctr){
    dv_nbrctr[dv]++;
    dl_at_dv[dv]=realloc(dl_at_dv[dv],sizeof(int)*dv_nbrctr[dv]);
}
void create_graph(){
    // dual graphs and links have prefix d
    // spatial and spatial-dual are s and sd
    int site;
    int vctr=0;
    lctr=0;
    dlctr=0;
    dvctr=0;
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
                wx_mark[lctr]=wx_smark[sl];
                wy_mark[lctr]=wy_smark[sl];
                if(i<3){
                    v_at_l[0][lctr] = v0;
                    v_at_l[1][lctr] = v1;
                }
                else {
                    v_at_l[0][lctr] = v0;
                    v_at_l[1][lctr] = v1;
                }
                add_l_to_v(v0,lctr);
                add_l_to_v(v1,lctr);
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
    //void stitch


}
int edge_notin_nbrs(int dv, int dl){
    int dv1=get_partner(dv,dl);
    int dvnbr;
    int rflag=1;
    for(int j=0; j<dv_nbrctr[dv]; j++){
        dvnbr= dv_nbr[j][dv1];
        if(dvnbr==dv1){
            rflag=0;
            break;
        }
    }
    return rflag;
}
void stich_in_time(){
    int sv, sdv, last, first;
    int dl;
    for(sv=0; sv<nsites; sv++){
        last=v_at_sv[sv];
        first=first_v[sv];
        if(first!=last){
            for(int j=0; j<dv_nbrctr[last]; j++){
                dl = dl_at_dv[j][sv];
                if (edge_notin_nbrs(first,dl)){


                }
            }
        }
    }
}






void remove_hyper_edges(){
    int l,dl,dl2,dv,l2;
    int maxj;
    int j,j2;
    for ( dv=0; dv<dvctr; dv++){
        j=0;
        while(j<dv_nbrctr[dv]){
            dl= dl_at_dv[dv][j];
            l=l_at_dl[dl];
            if (if_hyperedge[l]>1){
                j2=j+1;
                while(j2<dv_nbrctr[dv]){
                    dl2= dl_at_dv[dv][j];
                    l2=l_at_dl[dl2];
                    if (if_hyperedge[l2]==1){
                        dl_at_dv[dv][j]=l2;
                        dl_at_dv[dv][j2]=l;
                        dv_nbrctr[dv]--;
                        break;
                    }
                    j2++;
                }
                j=j2;
            }
            j++;
        }
    }
}
void add_dual_link( int sdv0, int sdv1){
    int sdl = get_sdl_from_sdv(sdv0, sdv1);
    int sl= get_sl_at_sdl(sdl);
    int l=l_at_sl[sl];
    //check if it corresponds to a hyperedge
    //if(if_hyperedge[l]==0){
    if(1){
        int dv0=dv_at_sdv[sdv0];
        int dv1=dv_at_sdv[sdv1];
        dv_at_dl[0][dlctr] = dv0;
        dv_at_dl[1][dlctr] = dv1;
        add_dl_to_dv(dv0,dlctr);
        add_dl_to_dv(dv1,dlctr);
        l_at_dl[dlctr]=l;
        dl_at_l[l]= dlctr;
        dlctr++;
        if_hyperedge[l]++;
    }
}


