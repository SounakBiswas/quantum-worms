#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"utils_dual_graph.h"
//graph data structures
vert **firstv ;
dvert **firstdv ;
vert **v_at_sv;
dvert **dv_at_sdv;
link **l_at_sl;

int lctr;
int dlctr;
int dvctr;
int vctr;
vert *get_partnerv(vert *v1,link *link);
dvert *get_partnerdv(dvert *dv1,dlink *dlink);

link *links;
dlink *dlinks;
vert *verts;
dvert *dverts;
void create_graph(){
  // dual graphs and links have prefix d
  // spatial and spatial-dual are s and sd
  int site;
  vctr=0;
  lctr=0;
  dlctr=0;
  dvctr=0;
  int op;
  int s1,s2,s3;
  int sdv0,sdv1,sdv2,sdv3;
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
        sdv0=dnbr[sdv][0];
        add_dual_link(sdv,sdv0);
        sdv1=dnbr[sdv][1];
        add_dual_link(sdv,sdv1);
        sdv2=dnbr[sdv][2];
        add_dual_link(sdv,sdv2);
      }

    }
  }

  for (int op_pos=0; op_pos <opstr_l; op_pos++) {
    op = opstr[op_pos];
    //new segment
    if(op < ntriangles){
      int sdv=op;
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

void stich_in_time(){
  int sv, sdv; 
  dvert *fdv,*ldv,*dv,*dv0,*dv1;
  vert *fv,*lv,*v0,*v1;
  dlink *dl,*dl0;
  link *l;
  int ix;
  for(sv=0; sv<nsites; sv++){
    lv= v_at_sv[sv];
    fv= firstv[sv];
    if(fv!=lv){
      for(int j=0; j<lv->nnbr; j++){
        l = lv->l[j];
        ix=get_link_ix(fv,l);
        if (ix==-1){
          v1= get_partnerv(lv,l);
          l->v0=fv;
          l->v1=v1;
          add_l_to_v(fv,l);
        }
        else
          remove_l_byix_from_v(v1,ix);

      }
    }
    remove_v(fv);
  }
  for(sdv=0; sdv<nsites; sdv++){
    ldv= dv_at_sdv[sdv];
    fdv= firstdv[sdv];
    if(fdv!=ldv){
      for(int j=0; j<ldv->nnbr; j++){
        dl = ldv->dl[j];
        ix=get_dlink_ix(fdv,dl);
        if (ix==-1){
          dv1= get_partnerdv(ldv,dl);
          dl->dv0=fdv;
          dl->dv1=dv1;
          add_dl_to_dv(fdv,dl);
        }
        else{
          dl0=fdv->dl[ix];
          (dl0->l)->he += (dl->l)->he;
          remove_dl_byix_from_dv(dv1,ix);
        }

      }
    }
    remove_dv(fdv);
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
    j2=dv->nnbr-1;
    while(j<dv->nnbr && j2>j){
      dl= dv->dl[j];
      l=dl->l;
      if (l->he >1){
        while(j2>j){
          dl2= dv->dl[j2];
          l2=dl2->l;
          if (l2->he==1){
            dv->dl[j]=dl2;
            dv->dl[j2]=dl;
            dv->nnbr--;
            break;
          }
          j2--;
        }
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


