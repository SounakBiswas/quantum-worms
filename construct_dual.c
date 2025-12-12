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
vert **v_at_sv;
dvert **dv_at_sdv;
int *sdv_at_dv;
link **l_at_sl;
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

link *links;
dlink *dlinks;
vert *verts;
dvert *dverts;

void add_l_to_v(vert *v, link* l){
    v->l=(link**)realloc(v->l,(v->nnbr+1)*sizeof(link *));
    v->l[v->nnbr]= l;
    v->nnbr++;
}
void add_dl_to_dv(dvert *dv, dlink* dl){
    dv->dl=(dlink**)realloc(dv->dl,(dv->nnbr+1)*sizeof(dlink *));
    dv->dl[dv->nnbr]= dl;
    dv->nnbr++;

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
    vert* v0,*v1;
    vert *v; link *l; 
    dvert *dv; dlink *dl; 
    dvert *dvp;
    int sv,sl;
    //add direct graph
    for (site=0; site<nsites; site++){
        v=&verts[vctr];
        v->s = sigma[site];
        v_at_sv[site]=v;
        vctr++;
    }
    int sv0,sv1;
    for (sv0=0; sv0<nsites; sv0++){
        for(int i=0; i<3; i++){
            l = & (links[lctr]);
            sv1 = nbr[site][i];
            sl = get_sl_from_sv(sv0,sv1);
            l_at_sl[sl] = l;
            l->wx=wx_smark[sl];
            l->wy=wy_smark[sl];
            l->v0=v_at_sv[sv0];
            l->v1=v_at_sv[sv1];
            l->d=(v0->s == v1->s)?1:0;
            l->he=0;
            lctr++;
        }
    }
    // add dual graph
    for (sdv=0; sdv<ndsites; sdv++){
        //sdv_at_dv[dvctr]=sdv;
        dv=&(dverts[dvctr]);
        dv_at_sdv[dsite]=dv;
        dv->ct=0;
        dvctr++;
    }
    for (sdv=0; dsite<ndsites; sdv++){

        for (int i=0; i<nplaqspersite; i++){
            //add dual links for uptriangles
            if(sdv%2==0){
                int sdv0=dnbr[sdv][0];
                add_dual_link(sdv,sdv0);
                int sdv1=dnbr[sdv][1];
                add_dual_link(sdv,sdv1);
                int sdv2=dnbr[sdv][2];
                add_dual_link(sdv,sdv2);
            }

        }
    }

    for (int op_pos=0; op_pos <opstr_l; op_pos++) {
        op = opstr[op_pos];
        //new segment
        if(op < ntriangles){
            int splaq=op;
            int sdv=splaq;
            dv = dv_at_sdv[sdv];
            dv->ct=1;
        }
        if (op >= ntriangles){
            sv0= (op-ntriangles)%nsites;
            if(op>ndiagops){
                sigma[sv0]*=-1;
            }

            v=&verts[vctr];
            v_at_sv[sv]=v;
            v->s=sigma[sv];
            v0=v_at_sv[sv];
            vctr++;
            //add links for new vertex (crucial to check for hyperedges
            for(int i=0; i<6; i++){
                l=&links[lctr];

                int sv1 = nbr[sv][i];
                sl=get_sl_from_sv(sv0,sv1);
                v1=v_at_sv[sv1];

                l_at_sl[sl] = l;
                l->wx=wx_smark[sl];
                l->wy=wy_smark[sl];
                l->v0=v0;
                l->v1=v1;
                l->d= (v0->s == v1->s) ? 1 : 0;
                l->he=0;
                add_l_to_v(v0,l);
                add_l_to_v(v1,l);
                lctr++;
            }
            //dual graph vertices
            //add 6 vertices
            for (int i=0; i<nplaqspersite; i++){
                dv= &dverts[dvctr];
                sdv=sdv_at_sv[sv][i];
                dv_at_sdv[sdv]=dv;
                dvctr++;
            }
            for (int i=0; i<nplaqspersite; i++){
                sdv=sdv_at_sv[sv][i];
                //uptriangles
                if(sdv%2==0){
                    int sdv0=dnbr[sdv][0];
                    add_dual_link(sdv,sdv0);
                    int sdv1=dnbr[sdv][1];
                    add_dual_link(sdv,sdv1);
                    int sdv2=dnbr[sdv][2];
                    add_dual_link(sdv,sdv2);
                }
                else{
                    int sdv2=dnbr[sdv][2];
                    add_dual_link(sdv2,sdv);
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
    link *l, *l2;
    int maxj;
    int j,j2;
    int dvix;
    dvert *dv;
    dlink * dl, *dl2;
    for ( dvix=0; dvix<dvctr; dvix++){
        dv= & dverts[dvix];
        j=0;
        while(j<dv->nnbr){
            dl= dv->dl[j];
            l=dl->l;
            if (l->he >1){
                j2=j+1;
                while(j2<dv->nnbr){
                    dl2= dv->dl[j2];
                    l2=dl2->l;
                    if (l2->he==1){
                        dv->dl[j]=dl2;
                        dv->dl[j2]=dl;
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
    link *l=l_at_sl[sl];
    dvert* dv0, *dv1;
    dlink *dl;
    //check if it corresponds to a hyperedge
    //if(if_hyperedge[l]==0){
    if(1){
        dl = &dlinks[dlctr];
        dv0=dv_at_sdv[sdv0];
        dv1=dv_at_sdv[sdv1];
        dl->dv0=dv0;
        dl->dv1=dv1;
        add_dl_to_dv(dv0,dl);
        add_dl_to_dv(dv1,dl);
        dl->l=l;
        l->he++;
        dlctr++;
    }
}


