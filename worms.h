int get_sdl_from_sdv(int sdv, int sdvnbr);
int get_sl_at_sdl(int sdl);
int get_sdl_at_sl(int sl);
int get_sl_from_sv(int sv,int svnbr);
void add_dual_link( int sdv0, int sdv1,int *dlctr);
extern int *sv_at_v;
extern int *v_at_sv;
extern int *dv_at_sdv;
extern int *sdv_at_dv;
extern int *l_at_sl;
extern int *dl_at_sdl;
extern int *l_at_dl;
extern int *dl_at_l;
extern int *p_at_sp;
extern int *v_at_l[2];
extern int *dv_at_dl[2];
extern int *dv2v[3];
extern int *dv_set;
extern int *dv_nbrctr;
extern int **dv_nbr; //nbr table for dual vertices
extern int **dl_at_dv; //site to link for dual vertices //int **dl_at_dv;
extern int *fsp;
extern int *if_hyperedge;
extern int *dimer;
extern int *backup_dimer;
extern int *dv_typ;

extern int vctr;
extern int lctr;
extern int dlctr;
extern int dvctr;

extern int *wx_mark;
extern int *wy_mark;
