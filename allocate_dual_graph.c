#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"worms.h"

void init_dual_graph(){
  int max_ndv=ndsites+6*(n_niop-n_triagop);
  int max_ndl=max_ndv*6;
  int max_nl=max_ndl;
  int max_nv=nsites+(n_niop-n_triagop);

  v_at_sv=(int*)malloc(nsites*sizeof(int));
  dv_at_sdv=(int*)malloc(ndsites*sizeof(int));
  l_at_sl = (int*) malloc(nbonds*sizeof(int));
  dl_at_sdl = (int*) malloc(nsdl*sizeof(int));

  sv_at_v=(int*)malloc(max_nv*sizeof(int));
  sdv_at_dv=(int*)malloc(max_ndv*sizeof(int));
  l_at_dl = (int*) malloc(max_ndl*sizeof(int));
  dl_at_l = (int*) malloc(max_nl*sizeof(int));
  for (int i =0; i<2; i++){
    v_at_l[i]=(int*)malloc(max_nl*sizeof(int));
    dv_at_dl[i]=(int*)malloc(max_ndl*sizeof(int));
  }
  for (int i =0; i<3; i++){
    dl_at_dv[i]=(int*)malloc(max_nv*sizeof(int));
  }
  fsp= (int*)malloc(max_nv*sizeof(int));
  if_hyperedge=(int*)malloc(max_nl*sizeof(int));
  dimer=(int*)malloc(max_nl*sizeof(int));;
  wx_mark = (int*)malloc(max_nl*sizeof(int));
  wy_mark = (int*)malloc(max_nl*sizeof(int));



  dv_typ=(int*)malloc(max_ndv*sizeof(int));

  l_at_v=(int**)malloc(max_nv*sizeof(int*));
  v_nbrctr = (int*)malloc(max_nv*sizeof(int));

  for(int i=0; i<max_nv; i ++){
    v_nbrctr[i]=6;
    l_at_v[i]=(int*)malloc(6*sizeof(int));
  }


}

void free_dual_graph(){
  int max_nv=nsites+(n_niop-n_triagop);

  free(v_at_sv);
  free(dv_at_sdv);
  free(l_at_sl );
  free(dl_at_sdl );

  free(sv_at_v);
  free(sdv_at_dv);
  free(l_at_dl );
  free(dl_at_l );
  for (int i =0; i<2; i++){
    free(v_at_l[i]);
    free(dv_at_dl[i]);
  }
  for (int i =0; i<3; i++){
    free(dl_at_dv[i]);
  }
  free(fsp);
  free(if_hyperedge);
  free(dimer);
  free(wx_mark );
  free(wy_mark );



  free(dv_typ);

  for(int i=0; i<max_nv; i ++){
    free(l_at_v[i]);
  }
  free(v_nbrctr );
  free(l_at_v);


}
