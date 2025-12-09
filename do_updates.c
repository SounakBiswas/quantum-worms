#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"global.h"
#include"mt19937ar.h"
#define max(a,b) ((a>b)?a:b)
void do_diag_withadjust(){
  double genrand_real2(void);
  int i;
  int site0,site1,site2;
  double k;
  int ransite,rantrngle;
  double P_add,P_del,P_site;
  k=beta*(nsites*h+ntriangles*2*jising);
  P_site=((nsites*beta*h)/k);
  for(i=0;i<opstr_l;i++) {
    if(opstr[i]==-1){//adding diagonal operators at identities
     P_add=k/(double)(opstr_l-n_niop);
      if(genrand_real2()<P_add){
	if(genrand_real2()< P_site){ //add a single site operator
        ransite=(int)(genrand_real2()*NSITES);
	  opstr[i]=ntriangles+ransite;
	  n_niop++;
	} 
	else{// attempt to add a triangle operator
	rantrngle=(int)(genrand_real2()*ntriangles);
	  site0=triangle[rantrngle][0];
	  site1=triangle[rantrngle][1];
	  site2=triangle[rantrngle][2];
	  if((sigma[site0]+sigma[site1]+sigma[site2])== 1 || (sigma[site0]+sigma[site1]+sigma[site2])== -1){// check if a triangle operator is allowed
	    opstr[i]=rantrngle;
	    n_niop++;
	    n_triagop++;
	  }
	}

      }
    }
    else if(opstr[i]<(ntrianglesandsites)){ 
      P_del=(double)(opstr_l-n_niop+1.0)/k;
      if(genrand_real2()<P_del){
	n_niop--;
	if(opstr[i]<(ntriangles)) n_triagop--;	
	opstr[i]=-1;
      }
    }  
    else  { //update spin state
      sigma[opstr[i]-ntrianglesandsites]*=-1;
    }     
  }
/* only exists in warmup */
  if(n_niop>(0.6*opstr_l)){// adjust length of operator string, should be done within the equilibration.
    new_opstr_l=opstr_l*1.5;
    opstr=(int *)realloc(opstr,new_opstr_l*sizeof(int));// reallocate memory
    divider=(int *)realloc(divider,new_opstr_l*sizeof(int));// reallocate memory
    for(i=opstr_l;i<new_opstr_l;i++)
      opstr[i]=divider[i]=-1;
    opstr_l=new_opstr_l;
  }
}//end of diagonal with adjust

void do_diag(){
  double genrand_real2(void);
  int i;
  int site0,site1,site2;
  double k;
  int ransite,rantrngle;
  double P_add,P_del,P_site;
  k=beta*(nsites*h+ntriangles*2*jising);
  P_site=((nsites*beta*h)/k);
  for(i=0;i<opstr_l;i++) {
    if(opstr[i]==-1){//adding diagonal operators at identities
     P_add=k/(double)(opstr_l-n_niop);
      if(genrand_real2()<P_add){
	if(genrand_real2()< P_site){ //add a single site operator
        ransite=(int)(genrand_real2()*NSITES);
	  opstr[i]=ntriangles+ransite;
	  n_niop++;
	} 
	else{// attempt to add a triangle operator
	rantrngle=(int)(genrand_real2()*ntriangles);
	  site0=triangle[rantrngle][0];
	  site1=triangle[rantrngle][1];
	  site2=triangle[rantrngle][2];
	  if((sigma[site0]+sigma[site1]+sigma[site2])== 1 || (sigma[site0]+sigma[site1]+sigma[site2])== -1){// check if a triangle operator is allowed
	    opstr[i]=rantrngle;
	    n_niop++;
	    n_triagop++;
	  }
	}

      }
    }
    else if(opstr[i]<(ntrianglesandsites)){ 
      P_del=(double)(opstr_l-n_niop+1.0)/k;
      if(genrand_real2()<P_del){
	n_niop--;
	if(opstr[i]<(ntriangles)) n_triagop--;	
	opstr[i]=-1;
      }
    }  
    else  { //update spin state
      sigma[opstr[i]-ntrianglesandsites]*=-1;
    }     
  }

}


int *label;
int labelmax=0;
int addlabel(){
  labelmax++;
  label[labelmax]=labelmax;
  return labelmax;

}
int find(int x){  // find representative element
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
int bind(int x, int y){
  int a=find(x);int b=find(y);
  if(a>b)  
    return label[a]=b;
  else
    return label[b]=a;

}

void do_clust(){//  The Hoshen Kopelman algorithm for Kedars algo
  int* flag;int *cluster;int *nlabels;int *clust_size;
  int i;int p;int j;
  int l,k;int c;
  int leg1,leg2,leg3;
  label=(int*)malloc((legs+1)*sizeof(int));
  cluster=(int*)malloc(legs*sizeof(int));
  nlabels=(int*)malloc((legs+1)*sizeof(int));
  clust_size=(int*)malloc((legs+1)*sizeof(int));
  double genrand_real2(void);
  labelmax=0;
  for(i=0;i<legs;i++){
    clust_size[i]=cluster[i]=nlabels[i]=0;
  }
  clust_size[legs]=nlabels[legs]=0;
  for(i=0;i<legs;i++){
      p=op_pos[i]; 
      if(opstr[p]>=2*nsites){
	if (cluster[link_list[i]]==0){
	  cluster[i]=addlabel(); // put in new cluster
	}
	else
	  cluster[i]=cluster[link_list[i]];
      }
      else{
	if(divider[p]!=-1){ // divided vertices
	  leg1=i+divider[p];
	  if (cluster[link_list[leg1]]==0){   
	    cluster[leg1]=cluster[leg1+3]=addlabel(); // put in new cluster
	  }
	  else
	    cluster[leg1]=cluster[leg1+3]=cluster[link_list[leg1]];
	  leg2=i+((divider[p]+1)%3);leg3=i+((divider[p]+2)%3);     
	  switch((!cluster[link_list[leg2]])+(!cluster[link_list[leg3]])){ // 4 leg part
	    case 2 :
	      cluster[leg2]=cluster[leg3]=cluster[leg2+3]=cluster[leg3+3]=addlabel();
	      break;
	    case 1 :
	      cluster[leg2]=cluster[leg3]=cluster[leg2+3]=cluster[leg3+3]=max(cluster[link_list[leg2]],cluster[link_list[leg3]]);
	      break;
	    case 0 :
	      cluster[leg2]=cluster[leg3]=cluster[leg2+3]=cluster[leg3+3]=bind(cluster[link_list[leg2]],cluster[link_list[leg3]]);
	      break;

	  }
	} 

	else{ // for the 6 leg vertices
	  switch((!cluster[link_list[i]])+(!cluster[link_list[i+1]])+(!cluster[link_list[i+2]])){
	    case 3 :
	      cluster[i]=cluster[i+1]=cluster[i+2]=cluster[i+3]=cluster[i+4]=cluster[i+5]=addlabel();
	      break;
	    case 2 :
	      if (cluster[link_list[i]]!=0)
		c=cluster[link_list[i]];
	      else if(cluster[link_list[i+1]]!=0)
		c=cluster[link_list[i+1]];
	      else
		c=cluster[link_list[i+2]];
	      cluster[i]=cluster[i+1]=cluster[i+2]=cluster[i+3]=cluster[i+4]=cluster[i+5]=c;
	      break;

	    case 1 :
	      if (cluster[link_list[i]]==0)
		c=bind(cluster[link_list[i+1]],cluster[link_list[i+2]]);
	      else if(cluster[link_list[i+1]]==0)
		c=bind(cluster[link_list[i]],cluster[link_list[i+2]]);
	      else
		c=bind(cluster[link_list[i+1]],cluster[link_list[i]]);
	      cluster[i]=cluster[i+1]=cluster[i+2]=cluster[i+3]=cluster[i+4]=cluster[i+5]=c;
	      break;
	    case 0 :
	      cluster[i+1]=cluster[i+2]= bind(cluster[link_list[i+1]],cluster[link_list[i+2]]);
	      cluster[i]=cluster[i+1]=cluster[i+2]=cluster[i+3]=cluster[i+4]=cluster[i+5]=bind(cluster[link_list[i]],cluster[i+1]);
	      break;

	  }
	}
	
	i=i+5;

      }
  
  }
  for(i=0;i<nsites;i++){ 
    if(first[i]!=-1) {
      cluster[first[i]]=cluster[last[i]]=bind(cluster[first[i]],cluster[last[i]]);
    }
  }
  labelmax=0;
  for(i=0;i<legs;i++){
      cluster[i]=find(cluster[i]);
      if(nlabels[cluster[i]]==0){
	labelmax++;nlabels[cluster[i]]=labelmax;
        cluster[i]=labelmax;
	clust_size[labelmax]++;
      }
      else{
	cluster[i]=nlabels[cluster[i]];
        clust_size[cluster[i]]++;
      }
  }
  if (take_stat_bool==1){
   int m;
   num_cluster+=labelmax;
  for(i=1;i<=labelmax;i++){
    m=(int)(N_CLUST_HIST_BINS*(((double)clust_size[i])/((double)legs)));
//    printf("%f\t%d\n",((double)clust_size[i]/(double)legs),m);
 //   getchar();
    hist_clust_size[m]++;
  }
  }
  free(nlabels);
  flag=(int*)malloc((labelmax+1)*sizeof(int));
  flag[0]=-2;
  for(i=1;i<=labelmax;i++)
    flag[i]=(genrand_real2()<0.5);
  for(i=0;i<legs;i++)
    if(flag[cluster[i]]==1){
      p=op_pos[i];
      // putchar('a');
      if((opstr[p]>=2*nsites)&&(opstr[p]<3*nsites))
	opstr[p]+=nsites;
      else if(opstr[p]>=3*nsites)
	opstr[p]-=nsites;

    }
  for(i=0;i<nsites;i++){
    if (first[i]==-1)
      sigma[i]=(genrand_real2()>0.5)?sigma[i]:-sigma[i];
    else 
      sigma[i]=(flag[cluster[first[i]]]==1)?-sigma[i]:sigma[i];
  }

  free(op_pos);
  free(link_list);
  free(cluster);
  free(flag);
  free(label);
  free(clust_size);
}



