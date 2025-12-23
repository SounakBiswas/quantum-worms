#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"fat.h"
//graph data structures
void stich_in_time();
void remove_hyper_edges();
int dvix;
vert **firstv ;
link **firstl ;
plaq **firstp ;
int lastsiteop[NSITES];
int dv_at_sdv[NDSITES];
vert **v_at_sv;
link **l_at_sl;
plaq **p_at_sdv;


vector links, dlinks, verts, dverts;
vector plaqs;
void create_graph(){
    // dual graphs and links have prefix d
    // spatial and spatial-dual are s and sd
    for(int i=0; i<NSITES; i++)
        lastsiteop[i]=-1;
    int site;
    int op;
    int s1,s2,s3;
    int sdv0,sdv1,sdv2,sdv3;
    int sdv;
    int sv0,sv1,sv2;
    vert* v0,*v1, *v2, *vold;
    vert *v; link *l; 
    int sl;
    dvix=0;
    //add direct graph
    for (sv0=0; sv0<nsites; sv0++){
        v=(vert*)top_ptr(&verts);
        v->id=(verts.top-1);
        v->s = sigma[sv0];
        v_at_sv[sv0]=v;
        firstv[sv0]=v;
        v->o1=-1;
        v->o1=-1;
    }
    for (sv0=0; sv0<nsites; sv0++){
        for(int i=0; i<3; i++){
            l= (link*)top_ptr(&links);
            l->id=links.top-1;
            sv1 = nbr[sv0][i];
            sl = get_sl_from_sv(sv0,sv1);
            firstl[sl]=l;
            l_at_sl[sl] = l;
            l->v0=v_at_sv[sv0];
            l->v1=v_at_sv[sv1];
            l->d=(l->v0->s == l->v1->s)?1:0;
            l->np=0;
        }
    }
    for (sdv=0; sdv<ndsites; sdv++){
        dv_at_sdv[sdv]=dvix;
        dvix++;
    }

    plaq *p,*lp;
    int sl0,sl1,sl2;
    int s0;
    for (int op_pos=0; op_pos <opstr_l; op_pos++) {
        op = opstr[op_pos];
        //new segment
        if(op!=-1 && op < ntriangles){
            sdv=op;
            lp=p_at_sdv[sdv];


            if(lp==NULL || lp->dvix!=dv_at_sdv[op]){
                int sdv=op;
                p=top_ptr(&plaqs);
                p->id=plaqs.top-1;
                //printf("adding, id=%d\n",p->id);
                p->dvix=dvix;
                //DEBUG_PRINT("sdv=%d\n",sdv);
                int i;
                for(i =0; i<3; i++){
                    sv0=triangle[op][i];
                    v0=v_at_sv[sv0];
                    p->v[i]=v0;
                    l=l_at_sl[ptolink[op][i]];
                    l->np++;
                    p->l[i]=l;
                    //printf("v %d\n",p->v[i]->id);

                }
                //printf("\n");
                p_at_sdv[sdv]=p;
                if(firstp[sdv]==NULL){
                    firstp[sdv]=p;
                }

            }

        }
        if (op!=-1 && op >= ntriangles){
            sv0= (op-ntriangles)%nsites;
            if(op>=ndiagops){
                sigma[sv0]*=-1;
            }
            vold=v_at_sv[sv0];
            vold->o2=op_pos;

            v= (vert*)top_ptr(&verts);
            v->id=verts.top-1;
            v_at_sv[sv0]=v;
            v->s=sigma[sv0];
            v->o1=op;
            v0=v_at_sv[sv0];
            lastsiteop[sv0]=op_pos;
            //add links for new vertex (crucial to check for hyperedges
            for(int i=0; i<6; i++){
                l= (link*) top_ptr(&links);
                l->id=links.top-1;

                int sv1 = nbr[sv0][i];
                sl=get_sl_from_sv(sv0,sv1);
                v1=v_at_sv[sv1];

                l_at_sl[sl] = l;
                l->v0=v0;
                l->v1=v1;
                l->d= (v0->s == v1->s) ? 1 : 0;
                l->np=0;
                //printf("added sl %d id %d \n",sl,l->id);
            }
            for (int i=0; i<nplaqspersite; i++){
                sdv=sdv_at_sv[sv0][i];
                dv_at_sdv[sdv]=dvix;
                dvix++;
                //DEBUG_PRINT("sdv=%d, ix=%d \n", sdv, dverts.top);
            }
        }
    }
    //void stitch
    //vert *fv,*lv;
    //int ix,ix1;
    plaq *p1,*p2;
    for(sdv=0; sdv<ndsites; sdv++){
        p=p_at_sdv[sdv];
        if(p!=NULL){
            assert(firstp[sdv]!=NULL);
            if(equal_plaqs(p,firstp[sdv])){
                p->id=-1;
            }
            else{
                for(int i=0; i<3; i++){
                    sv0=triangle[sdv][i];
                    v0=firstv[sv0];
                    p->v[i]=v0;
                    l=firstl[ptolink[op][i]];
                    l->np+= p->l[i]->np;
                    p->l[i]=l;
                }
            }
        }

    }
    vert *fv, *lv;
    for(sv0=0; sv0<nsites; sv0++){
        fv= firstv[sv0];
        fv->o1=lastsiteop[sv0];
        if(fv!=lv){
            lv->id=-1;
        }
    }
    print_graph();

}
int *label;
int addlabel();
int labelmax;
int find(int x){  // find representative element
    int y=x;int z; 
    while(label[y]!=y)
        y=label[y];
    while(x!=label[x]){// coalesce different labels into one
        int z=label[x];
        label[x]=y;
        x=z;
    }
    return y;

}
int bind(int x, int y){
    int a=find(x);int b=find(y);
    if(a>b)  
        return label[a]=b;
    else
        return label[b]=a;

}
int addlabel(){
    labelmax++;
    label[labelmax]=labelmax;
    return labelmax;

}
void new_cluster_update(){
    int  *cluster;
    int *nlabels;
    int *clust_size;
    label=(int*)malloc((legs+1)*sizeof(int));
    cluster=(int*)malloc(legs*sizeof(int));
    nlabels=(int*)malloc((legs+1)*sizeof(int));
    clust_size=(int*)malloc((legs+1)*sizeof(int));
    int pix;
    plaq *p;
    labelmax=0;
    vert *v0,*v1,*v2;
    vert *maj1,*maj2,*min;
    int ix0,ix1,ix2;;
    int ix;
    int id0,id1,id2;
    for(ix=0; ix <verts.top; ix++){
        cluster[ix]=-1;
        nlabels[ix]=-1;
    }
    for (pix=0; pix<plaqs.top; pix++){
        p = (plaq*)ix_ptr(&plaqs,pix);
        if(p->id!=-1){
            v0=p->v[0];
            v1=p->v[1];
            v2=p->v[2];
            for(int ix0=0; ix0<3; ix0++){
                ix1 = (ix0+1)%3;
                ix2 = (ix0+2)%3;
                if(p->v[ix0]->s == p->v[ix1]->s) {
                    maj1=p->v[ix0]; maj2=p->v[ix1];
                    min= p->v[ix2];
                    assert(maj1->id!=-1);
                    assert(maj2->id!=-1);
                    assert(min->id!=-1);
                    break;


                }
            }
            //cluster rule
            if(genrand_real2()<0.5){
                id0=maj2->id;
                id1=maj1->id;
                id2=min->id;
            }
            else {
                id0=maj1->id;
                id1=maj2->id;
                id2=min->id;
            }
            if(cluster[id0]==-1) {
                cluster[id0]=addlabel();
            }
            if(cluster[id1]==-1){ 
                cluster[id1]=cluster[id0];
            }
            else{
                cluster[id1]=cluster[id2]=bind(cluster[id1],cluster[id2]);
            }
        }
    }
    labelmax=0;
    for(ix=0;ix<verts.top;ix++){
        v0=(vert*)ix_ptr(&verts,ix);
        id0=v0->id;
        if(id0!=-1){
            cluster[id0]=find(cluster[id0]);
            if(nlabels[cluster[id0]]==-1){
                labelmax++;
                nlabels[cluster[id0]]=labelmax;
                cluster[id0]=labelmax;
                //clust_sidze[labelmax]++;
            }
            else{
                cluster[id0]=nlabels[cluster[id0]];
            }

        }
    }
    for(ix=0; ix <verts.top; ix++){
    }

}



