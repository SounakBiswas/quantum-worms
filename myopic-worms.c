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
  if(dv_typ[v]){
  while (xl=(int)(genrand_real2()*dv_nbrctr[v]), xl==nl);
  return xl;
    
  }

}

void myopic_worm(){
  int xl,nl;
  int start=(int)(genrand_real2()*dvctr);
  nl=-1;
  if(dv_typ[start]==1){
    for( int i=0; i<dv_nbrctr[start]; i++){
      if (dimer[dl_at_dv[start][i]]) {
        nl=dl_at_dv[start][i];
        break;
      }
    }
  }
  else if(dv_typ[start]==0){
    nl=dl_at_dv[start][(int)(genrand_real2()*dv_nbrctr[start])];
  }

  int v=get_partner(start,nl);



  


}
