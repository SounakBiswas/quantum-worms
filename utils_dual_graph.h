#include"worms.h"
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
int get_link_ix(vert *v0, link *l0);
int get_dlink_ix(dvert *dv0, dlink *dl0);

