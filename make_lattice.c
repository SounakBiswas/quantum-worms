#include<stdio.h>
#include"global.h"
#define  next_nbr(ii,jj,xx,yy)  ((xx+ii+lx)%lx + ((yy+jj+ly)%ly)*lx)

void make_lattice(void)
{
  int i,j,k,l,site,site1,site2,site3,sx,y;
  int orientation;
  for(i=0;i<lx;i++)
  {
    for(j=0;j<ly;j++)
    {
      site=i+lx*j;
       nbr[site][0] = next_nbr(i,j,1,0);
      nbr[site][3] = next_nbr(i,j,-1,0);
      nbr[site][2] = next_nbr(i,j,0,1);
      nbr[site][5] = next_nbr(i,j,0,-1);
      nbr[site][1] = next_nbr(i,j,1,1);
      nbr[site][4] = next_nbr(i,j,-1,-1);

      sublat[site]=(i+j)%3;


    }
  }
  for(i=0;i<(2*NSITES);i++){ // make triangles, 2*NSITES triangle
    site=i/2;
    orientation=i%2;
    if(orientation==0)  site1=nbr[site][0];  // up-pointing triangle
    else  site1=nbr[site][2]; //down-pointing triangle
    site2=nbr[site][1];

    triangle[i][sublat[site]]=site;
    triangle[i][sublat[site1]]=site1;
    triangle[i][sublat[site2]]=site2;

    if(orientation==0){
       dnbr[i][0]=2*nbr[site][0]+1;
       dnbr[i][2]=2*site+1;
       dnbr[i][2]=2*nbr[site][5]+1;
    }
    else{
      dnbr[i][0]=2*site;
      dnbr[i][1]=2*nbr[site][3];
      dnbr[i][2]=2*nbr[site][1];

    }

  }

}
