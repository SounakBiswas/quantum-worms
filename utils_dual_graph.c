#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"worms.h"
//stack
//push stack
void* top_ptr(vector *s){
  assert(s->top<s->size);
  void *elem= (void*)((char*)s->arr + s->top *s->size);
  s->top++;
  return elem;
}
void* ix_ptr(vector *s, int ix){
  assert(ix<=s->top);
  void *elem= (void*)((char*)s->arr + ix *s->size);
  return elem;
}
void alloc_vector(vector *s, int maxs){
  s->arr=malloc(maxs*s->size);
}
void free_stack(vector *s){
  free(s->arr);
}
void init_vector(vector *s, int max, size_t size){
    s->top=0;
    s->max=max;
    s->size=size;
    s->arr = (void*)malloc(max*size);
}
void remove_l_byix_from_v(vert *v, int ix){
  int i;
  v->l[ix]=v->l[v->nnbr];
  v->nnbr--;
  v->l[v->nnbr]=NULL;

}
void remove_dl_byix_from_dv(dvert *dv, int ix){
  //for( i=0; i<dv->nnbr; i++)
  //  if(dv->dl[i]==dl){
  //    break;
  //  }
  dv->dl[ix]=dv->dl[dv->nnbr];
  dv->nnbr--;
  dv->dl[dv->nnbr]=NULL;

}

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
void remove_v(vert *v0){
  free(v0->l);
  v0->nnbr=0;
}
void remove_dv(dvert *dv0){
  free(dv0->dl);
  dv0->nnbr=0;
}

void init_spatial_markers(){
  v_at_sv=(vert**)malloc(nsites*sizeof(vert*));
  firstv=(vert**)malloc(nsites*sizeof(vert*));
  firstdv=(vert**)malloc(ndsites*sizeof(vert*));
  dv_at_sdv=(dvert**)malloc(ndsites*sizeof(dvert*));
  l_at_sl = (link**) malloc(nbonds*sizeof(link*));

}
void free_spatial_markers(){
  free(v_at_sv);
  free(dv_at_sdv);
  free(l_at_sl  );
}

void init_dual_graph(){
  int max_ndv=ndsites+6*(n_niop-n_triagop);
  int max_ndl=max_ndv*6;
  int max_nl=max_ndl;
  int max_nv=nsites+(n_niop-n_triagop);
  links = (link*)malloc(max_nl*sizeof(link));
  verts = (vert*)malloc(max_nv*sizeof(link));
  dverts = (dvert*)malloc(max_ndv*sizeof(link));
  dlinks = (dlink*)malloc(max_ndl*sizeof(link));


}

void free_dual_graph(){
    int i;
    for(i=0; i<vctr; i++){
        free(verts[i].l);
    } 
    for(i=0; i<dvctr; i++){
        free(dverts[i].dl);
    } 
    free(verts);
    free(dverts);
    free(dlinks);
    free(links);

}
int get_dlink_ix(dvert *dv0, dlink *dl0){
  dlink* dl;
  dvert *dv2;
  int ix=-1;
  for (int j=0; j<dv0->nnbr; j++){
    dl=dv0->dl[j];
    if(dl==dl0){
      ix=j;
      break;
    }

  }
  return ix;
}
int get_link_ix(vert *v0, link *l0){
  link* l;
  vert *v2;
  int ix=-1;
  for (int j=0; j<v0->nnbr; j++){
    l=v0->l[j];
    if(l==l0){
      ix=j;
      break;
    }

  }
  return ix;
}
