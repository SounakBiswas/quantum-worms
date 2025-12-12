typedef struct link link;
typedef struct dlink dlink;
typedef struct vert vert;
typedef struct dvert dvert;

struct link {
  int v0;
  int v1;
  int dl;
  int he;
  int dimer;
};
struct dlink {
  int dv0;
  int dv1;
  int l;
};
struct vert {
  int nnbrs;
  int *links;
};
struct dvert {
  int nnbrs;
  int *dlinks;
};
