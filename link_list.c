#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
void makelinks(int update){ // type of updat : A,B or C coded as 0 1 or 2
  int i;FILE *fp;
  legs=2*n_niop + 4*n_triagop; // total number of legs
  op_pos=(int*)malloc(legs*sizeof(int));
 // leg_pos=(int*)malloc(legs*sizeof(int));
  link_list=(int*)malloc(legs*sizeof(int));
  int n_leg=0;
  int site0,site1,site2; // variables denoting sites for triangle operators
  for(i=0;i<nsites;i++)
    first[i]=last[i]=-1;
  for(i=0;i<opstr_l;i++)
    divider[i]=-1;
  for(i=0;i<opstr_l;i++){
    if(opstr[i]>=0){
    //  leg_pos[n_leg]=0;leg_pos[n_leg+1]=1;leg_pos[n_leg+2]=2;leg_pos[n_leg+3]=3; leg_pos[n_leg+4]=4; leg_pos[n_leg+5]=5; // update arrays that store position of leg in whatever vertex it belongs to, ie 1,2 3, or 4.
      if((opstr[i])<(NTRIANGLES)){ // if Triangle operator
      op_pos[n_leg]=op_pos[n_leg+1]=op_pos[n_leg+2]=op_pos[n_leg+3]= op_pos[n_leg+4]= op_pos[n_leg+5]=i; // update arrays to store corresponding postion in operator sequence for each leg
	site0=triangle[opstr[i]][0];
	site1=triangle[opstr[i]][1];
	site2=triangle[opstr[i]][2];
        if (first[site0]==-1)
          first[site0]=n_leg;
        if(last[site0]!=-1){
          link_list[n_leg]=last[site0];
	  link_list[last[site0]]=n_leg;
	}
        if (first[site1]==-1)
          first[site1]=n_leg+1;
        if(last[site1]!=-1){
          link_list[n_leg+1]=last[site1];
	  link_list[last[site1]]=n_leg+1;
	}


        if(first[site2]==-1)
          first[site2]=n_leg+2;
        if(last[site2]!=-1){
          link_list[n_leg+2]=last[site2];link_list[last[site2]]=n_leg+2;
        }
        last[site0]=n_leg+3;
        last[site1]=n_leg+4;
        last[site2]=n_leg+5;

	switch(update){
	 case 0:
	   if(sigma[site0]*sigma[site1]!=sigma[site0]*sigma[site2])
	     divider[i]=0;
           break;
	 case 1:
	   if(sigma[site1]*sigma[site0]!=sigma[site1]*sigma[site2])
	     divider[i]=1;
           break;
	 case 2:
	   if(sigma[site2]*sigma[site1]!=sigma[site2]*sigma[site0])
	     divider[i]=2;
           break;	   
	}
        n_leg=n_leg+6;
      }
      else{  // for single spin operators
	op_pos[n_leg]=op_pos[n_leg+1]=i;
//	op_pos[n_leg+2]=op_pos[n_leg+3]=op_pos[n_leg+4]=op_pos[n_leg+5]=i;
        site0=opstr[i]%nsites;
	if(opstr[i]>=ntrianglesandsites)
	  sigma[opstr[i]-ntrianglesandsites]*=-1;
        if(first[site0]==-1)
          first[site0]=n_leg;
        if(last[site0]!=-1){
          link_list[n_leg]=last[site0];
          link_list[last[site0]]=n_leg;
        }
        last[site0]=n_leg+1;
	//link_list[n_leg+2]=link_list[n_leg+3]=link_list[n_leg+4]=link_list[n_leg+5]=-1;
	n_leg=n_leg+2;
      }

    }

  } 
  for(i=0;i<nsites;i++){    // connect links across the boundary.
    if(last[i]!=-1){ 
      link_list[first[i]]=last[i];
      link_list[last[i]]=first[i];
    }
  }
/*    fp=fopen("data","w");
  for(i=0;i<legs;i++)
    fprintf(fp,"%d\t%d\t%d\t%d\n",i,link_list[i],op_pos[i],opstr[op_pos[i]]);
  fclose(fp);
  getchar();*/

}

