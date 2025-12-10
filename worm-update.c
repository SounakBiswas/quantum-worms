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
int get_xl(int v, int nl){
  int xl,dl;
  // 1 dimers
  int dn= dimer[l_at_dl[nl]];
  xl=nl;
  if(dv_typ[v]==1){
    if(dn==1){
      xl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]-1))];
      if (xl==nl)
        xl=dl_at_dv[v][dv_nbrctr[v]-1];
    }
    else if (dn==0){
      for( int i=0; i<dv_nbrctr[v]; i++){
        if (dl=dl_at_dv[v][i], dimer[l_at_dl[dl]]) {
          xl=dl;
          break;
        }
      }
    }
  }
  // 3 or one dimers allowed
  else if(dv_typ[v]==0){
      xl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]-1))];
  }
  return xl;

}

void myopic_worm(){
  int xl,nl;
  int start=(int)(genrand_real2()*dvctr);
  nl=-1;
  int n,x;
  int dl;
  int occ;
  int ln,lx;
  int dn,dx,v;
  nl=-1;
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

  v=get_partner(start,nl);
  while(v!=start){
     xl = get_xl(v, nl);
     nl=xl;
     v=get_partner(v,xl);
     ln=l_at_dl[nl];
     lx=l_at_dl[xl];
     dimer[ln]*=-1;
     dimer[lx]*=-1;
     //**check winding
  }

}
