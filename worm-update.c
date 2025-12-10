#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include"worms.h"
//nl and xl and entry and exit
int wx,wy;
int get_vertex_myopic();
int get_partner(int v1,int link){
  if(dv_at_dl[0][link]==v1)
    return dv_at_dl[0][link];
  else 
    return dv_at_dl[1][link];

}
int get_xl(int v, int nl, int dvtyp){
  int xl;
  // 1 dimers
  if(dv_typ[v]==1){
     if(dv_nbrctr[v]==1) 
       return nl;
     else{
       xl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]-1))];
       if (xl==nl)
         xl=dl_at_dv[v][dv_nbrctr[v]-1];
     }
  }
  // 3 or one dimers allowed
  else if(dv_typ[v]==0){
     if(dv_nbrctr[v]==1) 
       return nl;
     else{
       xl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]-1))];
       if (xl==nl)
         xl=dl_at_dv[v][dv_nbrctr[v]-1];
     }
  }
  else if(dv_typ[v]==2){
     if(dimer[l_at_dl[nl]]==1) 
       return nl;
     else{
       //**assuming hyper edge connections have been removed
       xl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]))];
     }
  }
  return xl;

}

void myopic_worm(){
  int xl,nl;
  int start=(int)(genrand_real2()*dvctr);
  nl=-1;
  int n,x;
  int dl;
  if(dv_typ[start]==1){
    for( int i=0; i<dv_nbrctr[start]; i++){
      if (dl=dl_at_dv[start][i], dimer[l_at_dl[dl]]) {
        nl=dl;
        break;
      }
    }
  }
  else if(dv_typ[start]==0){
    nl=dl_at_dv[start][(int)(genrand_real2()*dv_nbrctr[start])];
  }

  int v=get_partner(start,nl);
  n=start;

  //do the worm pivot
  if(dv_typ[v]==1){
    for( int i=0; i<dv_nbrctr[v]; i++){
      if (dimer[dl_at_dv[start][i]]) {
        nl=dl_at_dv[start][i];
        break;
      }
    }
  }




  


}
