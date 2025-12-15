#include<stdio.h>
#include"global.h"
#include"dual_graphs.h"
#define  next_nbr(ii,jj,xx,yy)  ((xx+ii+lx)%lx + ((yy+jj+ly)%ly)*lx)
// 2     1
// --\  /
//3___\/__0
//    /\
//   /  \
//  4    5
int get_sl_from_sv(int sv1,int sv2){
    int i;
    int sl=-1;
    for(i=0;i<3;i++){
        if(nbr[sv1][i]==sv2){
            sl=3*sv1+i;
        }
    }
    if(sl!=-1)
        return sl;
    else{
        for(i=0;i<3;i++){
            if(nbr[sv2][i]==sv1){
                sl=3*sv2+i;
            }
        }
    }
    return sl;



}
int get_sdl_from_sdv(int sdv1,int sdv2){
    int i;
    int sl=-1;
    if(sdv1%2==0){
        for(i=0; i<3; i++){
            if (dnbr[sdv1][i]==sdv2){
                return 3*sdv1+i;
            }
        }
    }
    else{
        for(i=0; i<3; i++){
            if (dnbr[sdv2][i]==sdv1){
                return 3*sdv2+i;
            }
        }
    }
    return -1;
}

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

            sdv_at_sv[site][0] = 2*site;
            sdv_at_sv[site][1] = 2*site+1;
            sdv_at_sv[site][2] = 2*nbr[site][3];
            sdv_at_sv[site][3] = 2*nbr[site][4]+1;
            sdv_at_sv[site][4] = 2*nbr[site][4];
            sdv_at_sv[site][5] = 2*nbr[site][5]+1;


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

    }
    for(i=0; i<NDSITES; i++){
        site=i/2;
        orientation=i%2;
        if(orientation==0){
            dnbr[i][0]=2*nbr[site][0]+1;
            dnbr[i][1]=2*site+1;
            dnbr[i][2]=2*nbr[site][5]+1;
        }
        else{
            dnbr[i][1]=2*site;
            dnbr[i][0]=2*nbr[site][3];
            dnbr[i][2]=2*nbr[site][1];

        }
    }

    int dlix,upsite;
    int nuptriangles=nsites;
    for(i=0; i<nuptriangles; i++){
        sl_at_sdl[3*i]= 3*nbr[upsite][0]+2;
        sl_at_sdl[3*i+1]= 3*upsite+1;
        sl_at_sdl[3*i+2]= 3*upsite;
    }

}
