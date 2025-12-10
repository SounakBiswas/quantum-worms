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
int get_xdl(int v, int nl){
  int xdl,dl;
  // 1 dimers
  int dn= dimer[l_at_dl[nl]];
  xdl=nl;
  if(dv_typ[v]==1){
    if(dn==1){
      xdl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]-1))];
      if (xdl==nl)
        xdl=dl_at_dv[v][dv_nbrctr[v]-1];
    }
    else if (dn==0){
      for( int i=0; i<dv_nbrctr[v]; i++){
        if (dl=dl_at_dv[v][i], dimer[l_at_dl[dl]]) {
          xdl=dl;
          break;
        }
      }
    }
  }
  // 3 or one dimers allowed
  else if(dv_typ[v]==0){
      xdl=dl_at_dv[v][(int)(genrand_real2()*(dv_nbrctr[v]-1))];
  }
  return xdl;

}

void worm_update(){
  int xdl,ndl;
  int start=(int)(genrand_real2()*dvctr);
  ndl=-1;
  int n,x;
  int dl;
  int occ;
  int ln,lx;
  int dn,dx,v;
  ndl=-1;
  //counters to check winding number parities
  int wx_ctr=0;
  int wy_ctr=0;
  if(dv_typ[start]==1){
    for( int i=0; i<dv_nbrctr[start]; i++){
      if (dl=dl_at_dv[start][i], dimer[l_at_dl[dl]]) {
        ndl=dl;
        break;
      }
    }
  }
  else if(dv_typ[start]==0){
    ndl=dl_at_dv[start][(int)(genrand_real2()*dv_nbrctr[start])];
  }

  v=get_partner(start,ndl);
  ln=l_at_dl[ndl];
  wx_ctr+=wx_mark[ln];
  wy_ctr+=wy_mark[ln];
  while(v!=start){
     xdl = get_xdl(v, ndl);
     ndl=xdl;
     v=get_partner(v,xdl);
     ln=l_at_dl[ndl];
     lx=l_at_dl[xdl];
     dimer[ln]*=-1;
     dimer[lx]*=-1;
     //**check winding
     //**how to check winding ?
     wx_ctr+=wx_mark[lx];
     wy_ctr+=wy_mark[lx];
  }

}
