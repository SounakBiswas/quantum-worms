#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#include <assert.h>
void canonical(){
  long int i ;  
  int current_leg,start_leg;
  int pos1,pos2,x,y;
  int ctr,fail,count;

  int next_leg,first_leg;
  int * clus;
  clus=(int *)malloc((n_niop+1)*sizeof(int));
  int update_ctr;
  int op,op_p,bond;
  int site;
  int direction;
  int min_site;
  n_flipop=0;
  for(i=0; i<opstr_l; i++){
      if (opstr[i]>ntrianglesandsites)
          n_flipop++;
  }
  //printf("n_niop %d n_triagop %d flips %d wolffsteps=%d \n",n_niop,n_triagop,n_flipop,wolffsteps);
  int trian;

  if(n_niop!=n_triagop){

    for(update_ctr=0;update_ctr<wolffsteps;update_ctr++){
      i =(int)(genrand_real2()*2.0*(n_niop-n_triagop));
      start_leg= current_leg =pos_gammalegs[i];
      direction=start_leg%2?-1:1;
      site=(opstr[op_pos[start_leg]]-ntriangles)%nsites;
      op=opstr[op_pos[start_leg]];
      clus[0]=op_pos[start_leg];
      //printf("start pos op %d %d \n",op_pos[start_leg],opstr[clus[0]]);
      ctr=0;
      fail=0;
      count=0;
      //  printf("nsites %d \n",nsites);
      //make 1-D clus
      while(1){
        ctr++;
        next_leg=link_list[current_leg];
        op_p=op_pos[next_leg];
        clus[ctr]=op_p;
        if(opstr[op_p]>=ntriangles)  //loop terminates
          break;
        trian = opstr[op_p]; 
        min_site=triangle[trian][minority[op_p]];
        if(min_site==site){ //loop fails
          fail=1;
          break;
        }

        current_leg=next_leg-3*direction;
      }
      //printf("ctr=%d, start=%d, rand=%f, wt=%f, count=%d, fail=%d, accept=%d\n",ctr,start_leg,rand,wt,count,fail,(rand<wt)&&(fail==0));
      //printf("fail=%d\n",fail);
      //getchar();

      if(fail==0){
       //printf("accept flip\n");
        //update starting operator
        op=opstr[clus[0]];
        //printf("%d %d\n",op,clus[0]);
        if(op>=ntrianglesandsites)
          op-=nsites;
        else op+=nsites;
        opstr[clus[0]]=op;

        //update ending operator
        op=opstr[clus[ctr]];
        if(op>=ntrianglesandsites)
          op-=nsites;
        else 
          op+=nsites;
        opstr[clus[ctr]]=op;

        //update string of J operators
        for(i=1;i<ctr;i++){
          op_p=clus[i];
          op=opstr[op_p];
          minority[op_p]=3-(sublat[site]+minority[op_p]);

        }

        //flip spin, if applicable
        if((start_leg<next_leg)&&(direction==1))
          sigma[site]*=-1;
        if((start_leg>next_leg)&&(direction==-1))
          sigma[site]*=-1;

      }
    }
   }
  free(clus);


}

