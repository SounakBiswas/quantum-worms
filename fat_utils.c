#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include<assert.h>
#include"fat.h"
//stack
//push stack
//returns top pointer and advances
void* top_ptr(vector *s){
   assert((s->top) < (s->max));
  void *elem= (void*)((char*)(s->arr) + (s->top) *(s->size) );
  s->top++;
  return elem;
}
void* ix_ptr(vector *s, int ix){
  assert(ix<=s->top);
  void *elem= (void*)((char*)(s->arr) + ix *(s->size));
  return elem;
}
void alloc_vector(vector *s, int maxs){
  s->arr=malloc(maxs*s->size);
}
void free_vector(vector *s){
  free(s->arr);
}
void init_vector(vector *s, int max, size_t size){
    s->top=0;
    s->max=max;
    s->size=size;
    s->arr = (void*)malloc(max*size);
}

void init_spatial_markers(){
  v_at_sv=(vert**)malloc(nsites*sizeof(vert*));
  firstv=(vert**)malloc(nsites*sizeof(vert*));
  firstl=(link**)malloc(nbonds*sizeof(link*));
  firstp=(plaq**)malloc(ndsites*sizeof(plaq*));
  l_at_sl = (link**) malloc(nbonds*sizeof(link*));
  p_at_sdv = (plaq**) malloc(ndsites*sizeof(plaq*));
  for(int i=0; i<ndsites; i++){
      p_at_sdv[i]=NULL;
      firstp[i]=NULL;
  }

}
void free_spatial_markers(){
  free(v_at_sv);
  free(p_at_sdv);
  free(l_at_sl  );
  free(firstv);
  free(firstl);
  free(firstp);
}

void init_dual_graph(){
  int max_ndv=ndsites+6*(n_niop-n_triagop);
  int max_ndl=max_ndv*6;
  int max_nl=max_ndl;
  int max_nv=nsites+(n_niop-n_triagop);
  init_vector(&links, max_nl, sizeof(link));
  init_vector(&verts, max_nv, sizeof(vert));
  init_vector(&plaqs,max_ndv, sizeof(plaq));



}

void free_dual_graph(){

    free_vector(&links);
    free_vector(&verts);
    free_vector(&plaqs);

}
int equal_plaqs(plaq *p1, plaq *p2){
    for (int i =0; i<3; i++)
        if(p1->v[i]!=p2->v[i])
            return 0;
    return 1;

}
void print_graph(){
    int i;
    plaq *p;
    for(i=0; i<plaqs.top; i++){
        p=(plaq*)ix_ptr(&plaqs,i);
        printf("plaq, id=%d, dvix=%d\n",p->id, p->dvix);
        if(p->id!=-1){
            printf("verts %d %d %d \n",p->v[0]->id,p->v[1]->id,p->v[2]->id);
            printf("links %d(%d) %d(%d) %d(%d) \n",p->l[0]->id,p->l[0]->np,p->l[1]->id,p->l[1]->np,p->l[2]->id,p->l[2]->np);
        }
        //else{
        //    printf("removed plaq \n");
        //}


    }
}
