#include<stdlib.h>
#include<stdio.h>
int get_sdl_from_sdv(int sdv, int sdvnbr);
int get_sl_at_sdl(int sdl);
int get_sl_from_sv(int sv,int svnbr);
void print_dual();
void fprint_dual();
void add_dual_link( int sdv0, int sdv1);
extern int *sv_at_v;


extern int *sdv_at_dv;
extern int *dl_at_sdl;
extern int *l_at_dl;
extern int *dl_at_l;
extern int **dv_nbr; //nbr table for dual vertices
//extern int *dl_at_dv[3]; //site to link for dual vertices //int **dl_at_dv;



typedef struct link link;
typedef struct dlink dlink;
typedef struct vert vert;
typedef struct dvert dvert;
typedef struct plaq plaq;
typedef struct vector vector;
struct vector{
  void* arr;
  int top;
  int max;
  size_t size;
};

struct link{
  vert *v0, *v1;
  int wx,wy;
  dlink *dl;
  int d;
  int he;
  int id;
};
struct dlink{
  dvert *dv0, *dv1;
  link *l;
  int id;
};
struct vert{
  int nnbr;
  link **l;
  int s;
  int id;
};
struct dvert{
  int nnbr;
  dlink **dl;
  int ct;
  int id;
};
struct plaq{
  vert *v0,*v1,*v2;
  link *l0, *l1, *l2;
  int id;
};
extern vert **v_at_sv;
extern link **l_at_sl;
extern dvert **dv_at_sdv;
extern vert **firstv ;
extern dvert **firstdv ;

extern vector links,verts,dlinks,dverts;
extern int lctr;
extern int dlctr;
extern int vctr;
extern int dvctr;
void init_spatial_markers();
void free_spatial_markers();
void init_dual_graph();
void free_dual_graph();
void remove_l_byix_from_v(vert *v, int ix);
void remove_dl_byix_from_dv(dvert *dv, int ix);

void add_l_to_v(vert *v, link* l);
void add_dl_to_dv(dvert *dv, dlink* dl);
void remove_v(vert *v0);
void remove_dv(dvert *dv0);
int get_link_ix(vert *v0, vert *l0);
int get_dlink_ix(dvert *dv0, dvert *dl0);
void* top_ptr(vector *s);
void* ix_ptr(vector *s, int ix);
