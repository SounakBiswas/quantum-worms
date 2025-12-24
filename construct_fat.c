#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h" 
#include"fat.h"
#include"time.h"
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
    time_t t0=clock();
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
        v->o2=-1;
        v->sv=sv0;
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
                DEBUG_PRINT("adding, id=%d at sdv %d\n",p->id,sdv);
                if(lp!=NULL){
                    DEBUG_PRINT(" old dvix=%d (id=%d) at sdv %d\n",lp->dvix,lp->id,dv_at_sdv[op]);
                }
                else
                    DEBUG_PRINT("first here at %d \n",dv_at_sdv[op] );
                p->dvix=dv_at_sdv[op];
                //DEBUG_PRINT("sdv=%d\n",sdv);
                int i;
                for(i =0; i<3; i++){
                    sv0=triangle[op][i];
                    v0=v_at_sv[sv0];
                    p->v[i]=v0;
                    //add_p_to_v(v0,p);
                    assert(v0->id!=-1);
                    DEBUG_PRINT("v %d\n",p->v[i]->id);
                }
                //if(p->v[0]->s == p->v[1]->s){
                //    p->min=2;
                //}
                //else if(p->v[0]->s == p->v[2]->s){
                //    p->min=1;
                //}
                //else
                //    p->min=0;
                DEBUG_PRINT(" added sigma, %d %d %d \n",sigma[p->v[0]->sv],sigma[p->v[1]->sv],sigma[p->v[2]->sv]);
                DEBUG_PRINT(" sp at v %d %d %d, sum=%d \n",p->v[0]->s,p->v[1]->s,p->v[2]->s,p->v[0]->s + p->v[1]->s + p->v[2]->s);
                assert(abs(p->v[0]->s + p->v[1]->s + p->v[2]->s)==1);
                //DEBUG_PRINT("\n");
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
            v->o1=op_pos;
            v->o2=-1;
            v->sv=sv0;
            v0=v_at_sv[sv0];
            lastsiteop[sv0]=op_pos;
            //add links for new vertex (crucial to check for hyperedges
            for (int i=0; i<nplaqspersite; i++){
                sdv=sdv_at_sv[sv0][i];
                dv_at_sdv[sdv]=dvix;
                dvix++;
            }
        }
    }
    //void stitch
    //vert *fv,*lv;
    //int ix,ix1;
    time_t t1=clock();
    plaq *p1,*p2;
    int pix;
    vert *lv,*fv;
    for(pix=0; pix<plaqs.top; pix++){
        p= (plaq*)ix_ptr(&plaqs,pix);
        for(int i=0; i<3; i++){
            v0=p->v[i];
            sv0= p->v[i]->sv;
            lv=v_at_sv[sv0];
            fv=firstv[sv0];
            if(v0==lv && lv!=fv){
                p->v[i]=fv;
            }
        }

    }
    for(sv0=0; sv0<nsites; sv0++){
        fv= firstv[sv0];
        fv->o1=lastsiteop[sv0];
        lv=v_at_sv[sv0];
        DEBUG_PRINT("sv0=%d, fv=%d (%p), lv=%d (%p)\n",sv0,fv->id,fv,lv->id,lv);
        if(lv!=fv){
            lv->id=-1;
        }
    }
    tstitch+=(double)(clock()-t1)/(1.0*CLOCKS_PER_SEC);
    //print_graph();
    tmake+=(double)(clock()-t0)/(1.0*CLOCKS_PER_SEC);

}
int *label;
int labelmax;
int  *cluster;
static int find(int x){  // find representative element
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
static int addlabel(){
    labelmax++;
    label[labelmax]=labelmax;
    return labelmax;

}
static int bind(int x, int y){
    if(x==-1){
        x=addlabel();
    }
    if(y==-1){
        y=addlabel();
    }
    int a=find(x);int b=find(y);
    if(a>b)  
        return label[a]=b;
    else
        return label[b]=a;

}
void flip_at_op_pos(int p){
    if((opstr[p]>=2*nsites)&&(opstr[p]<3*nsites))
        opstr[p]+=nsites;
    else if(opstr[p]>=3*nsites)
        opstr[p]-=nsites;
}
void split_plaq_rand(plaq *p);
void split_plaq_prem(plaq *p,int r);

void new_cluster_update(int typ){
    time_t t0=clock();
    int *nlabels;
    int *clust_size;
    label=(int*)malloc((verts.top)*sizeof(int));
    cluster=(int*)malloc(verts.top*sizeof(int));
    nlabels=(int*)malloc((verts.top)*sizeof(int));
    clust_size=(int*)malloc((verts.top)*sizeof(int));
    int pix;
    plaq *p;
    labelmax=-1;
    vert *v0,*v1,*v2;
    vert *maj1,*maj2,*min;
    int ix0,ix1,ix2;;
    int ix;
    int id0,id1,id2;
    for(ix=0; ix <verts.top; ix++){
        v0=(vert*)ix_ptr(&verts,ix);
        if(v0->id !=-1)
          cluster[ix]=addlabel();
        else
          cluster[ix]=-1;

        nlabels[ix]=-1;

    }
    for (pix=0; pix<plaqs.top; pix++){
        p = (plaq*)ix_ptr(&plaqs,pix);
        if(typ==3){
           split_plaq_rand(p);
        }
        else{
            assert(typ<3);
           split_plaq_prem(p,typ);
        }
    }
    //*stitch in update*//
    //for(ix=0; ix<nsites;ix++){
    //    id0 = firstv[ix]->id;
    //    id1 = v_at_sv[ix]->id;
    //    if(id0!=id1){
    //        cluster[id0]=cluster[id1]=bind(cluster[id0],cluster[id1]);
    //    }
    //}
    //***/
    labelmax=-1;
    for(ix=0;ix<verts.top;ix++){
        v0=(vert*)ix_ptr(&verts,ix);
        id0=v0->id;
        if( id0!=-1 && cluster[id0]!=-1){
            cluster[id0]=find(cluster[id0]);
            if(nlabels[cluster[id0]]==-1){
                labelmax++;
                DEBUG_PRINT("lmax=%d, id0=%d\n",labelmax,id0);
                nlabels[cluster[id0]]=labelmax;
                cluster[id0]=labelmax;
                //clust_sidze[labelmax]++;
            }
            else{
                cluster[id0]=nlabels[cluster[id0]];
            }

        }
    }
    free(nlabels);
    //DEBUG_PRINT("num clust= %d \n",labelmax+1);
    //getchar();
    //cluster check
    for (pix=0; pix<plaqs.top; pix++){
        p = (plaq*)ix_ptr(&plaqs,pix);
        if(typ<3 && p->id!=-1){
        ix0=typ;
        ix1 = (ix0+1)%3;
        ix2 = (ix0+2)%3;
        v0=p->v[ix0];
        v1=p->v[ix1];
        v2=p->v[ix2];
        id0=p->v[ix0]->id;
        id1=p->v[ix1]->id;
        id2=p->v[ix2]->id;

        if( (v0->s == v1->s) || (v0->s == v2->s)){
            assert(cluster[id1]==cluster[id2]);
            //DEBUG_PRINT("ids(s): %d %d %d \n",id0,id1,id2);
            //DEBUG_PRINT("clu(s): %d %d %d \n",cluster[id0],cluster[id1],cluster[id2]);
        }
        
        else{
            assert(v1->s == v2->s);
            //DEBUG_PRINT("ids(f): %d %d %d \n",id0,id1,id2);
            //DEBUG_PRINT("clu(f): %d %d %d \n",cluster[id0],cluster[id1],cluster[id2]);
            assert(cluster[id0]== cluster[id1]);
            assert(cluster[id0]== cluster[id2]);
            assert(cluster[id2]== cluster[id1]);
            //DEBUG_PRINT("freeze: %d %d(%d) %d %d(%d) %d %d(%d) \n",id0,cluster[id0],label[cluster[id0]], id1,cluster[id1],label[cluster[id1]], id2,cluster[id2],label[cluster[id2]]);


        }
            
        }
    }
    int *flag=(int*)malloc((labelmax+1)*sizeof(int));
    for(int i=0;i<=labelmax;i++)
        flag[i]=(genrand_real2()<0.5);
    for(ix=0; ix <verts.top; ix++){
        v0=(vert *)ix_ptr(&verts,ix);
        id0=v0->id;
        if(id0 !=-1){
            if (cluster[id0]!=-1 && flag[cluster[id0]]==1){
                if ( v0->o1 !=-1){
                    assert(v0->o2!=-1);
                    flip_at_op_pos(v0->o1);
                    flip_at_op_pos(v0->o2);
                    //v0->s *= -1;
                }
                //if(v0->o1 != -1)
                //    flip_at_op_pos(v0->o1);
                //if(v0->o2 != -1)
                //    flip_at_op_pos(v0->o2);



            }

        }
    }
    vert *fv, *lv;
    for(ix=0; ix <nsites; ix++){
        fv=firstv[ix];
        id0=fv->id;
        assert(id0!=-1);
        if(cluster[id0]==-1){
            sigma[ix]=(genrand_real2()<0.5)?1:-1;
            DEBUG_PRINT("i=%d, clus=free\n",ix);
            //fv->s *=-1;
        }
        else{
            DEBUG_PRINT("i=%d, clus=%d\n",ix,cluster[id0]);
            if(flag[cluster[id0]]==1){
                sigma[ix]*=-1;
            }
        } 

    }
    //getchar();
    DEBUG_PRINT("******************************\n");
    free(cluster);
    free(flag);
    free(label);
    free(clust_size);
    DEBUG_PRINT("update compeleted \n");
    tclust+=(double)(clock()-t0)/(1.0*CLOCKS_PER_SEC);

}


void fat_update(int typ){
    init_spatial_markers();
    //DEBUG_PRINT("spatial markers done\n");
    init_dual_graph();
    //DEBUG_PRINT("dual allocations done\n");
    create_graph();
    new_cluster_update(typ);
    //DEBUG_PRINT("graph created\n");
    free_dual_graph();
    free_spatial_markers();
}

void split_plaq_rand(plaq* p){
    vert *v0, *v1, *v2;
    int ix0, ix1, ix2;
    int id0, id1, id2;
    vert *maj1,*maj2,*min;
    if(p->id!=-1){
        DEBUG_PRINT("pid=%d\n",p->id);
        DEBUG_PRINT("verts=%d %d %d\n",p->v[0]->id, p->v[1]->id,p->v[2]->id);
        v0=p->v[0];
        v1=p->v[1];
        v2=p->v[2];
        assert(v0->id!=-1);
        assert(v1->id!=-1);
        assert(v2->id!=-1);
        for(int ix0=0; ix0<3; ix0++){
            ix1 = (ix0+1)%3;
            ix2 = (ix0+2)%3;
            if(p->v[ix0]->s == p->v[ix1]->s) {
                maj1=p->v[ix0]; maj2=p->v[ix1];
                min= p->v[ix2];
                assert(maj1->id!=-1);
                assert(maj2->id!=-1);
                assert(min->id!=-1);
                DEBUG_PRINT(" clust build :%d %d %d \n",maj1->s, maj2->s, min->s);
                assert(abs(maj1->s + maj2->s + min->s)==1);
                assert(maj1->s * maj2->s ==1);
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
        cluster[id1]=cluster[id2]=bind(cluster[id1],cluster[id2]);
    }
}
void split_plaq_prem(plaq* p, int typ){
    vert *v0, *v1, *v2;
    int ix0, ix1, ix2;
    int id0, id1, id2;
    vert *maj1,*maj2,*min;
    if(p->id!=-1){
        DEBUG_PRINT("pid=%d\n",p->id);
        DEBUG_PRINT("verts=%d %d %d\n",p->v[0]->id, p->v[1]->id,p->v[2]->id);
        v0=p->v[0];
        v1=p->v[1];
        v2=p->v[2];
        assert(v0->id!=-1);
        assert(v1->id!=-1);
        assert(v2->id!=-1);
        ix0=typ;
        ix1 = (ix0+1)%3;
        ix2 = (ix0+2)%3;
        v0=p->v[ix0];
        v1=p->v[ix1];
        v2=p->v[ix2];

        if( (v0->s == v1->s) || (v0->s == v2->s)){
            id0=v0->id;
            id1=v1->id;
            id2=v2->id;
            if(cluster[id0]==-1) {
                cluster[id0]=addlabel();
            }
            cluster[id1]=cluster[id2]=bind(cluster[id1],cluster[id2]);
        }
        else{
            assert(v1->s == v2->s);
            id0=v0->id;
            id1=v1->id;
            id2=v2->id;
            cluster[id0]=cluster[id1]=bind(cluster[id0],cluster[id1]);
            cluster[id0]=cluster[id2]=cluster[id1]=bind(cluster[id1],cluster[id2]);
            //DEBUG_PRINT("freeze: %d %d(%d) %d %d(%d) %d %d(%d) \n",id0,cluster[id0],label[cluster[id0]], id1,cluster[id1],label[cluster[id1]], id2,cluster[id2],label[cluster[id2]]);


        }

        //cluster rule
    }
}
