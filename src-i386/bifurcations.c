#include <R.h>
#include <Rinternals.h>

// find the connecting edge from a node
int getEdg(int nd, int * edge1, int * edge2,
           int nedges)
{
  int i;
  int edg=-1;
  for (i=0;i<nedges;i++) {
    if (nd==edge2[i]) {
      edg=i;
      break;
    }
  }
  return edg;
}

// find the connecting edge from an edge
int cnnctngEdg(int edge, int * edge1, int * edge2,
               int nedges)
{
  int nxt;
  int nd=edge1[edge];
  nxt = getEdg(nd, edge1, edge2, nedges);
  return nxt;
}

// return array of edge bipartitions
SEXP bifurcations(SEXP e1, SEXP e2,
                  SEXP nts)
{
  int i;
  int edg;
  SEXP res;
  int* edge1 = INTEGER(e1);
  int* edge2 = INTEGER(e2);
  int nedges = length(e1);
  int ntips = asInteger(nts);
  PROTECT(res=allocMatrix(INTSXP, nedges, ntips));
  int n = length(res);
  for(i=0;i<n; i++) {
    INTEGER(res)[i] = 0;
  }
  for(i=0;i<ntips; i++) {
    edg=getEdg(i+1, edge1, edge2, nedges);
    INTEGER(res)[edg + i * nedges]=1;
    edg=cnnctngEdg(edg, edge1, edge2, nedges);
    while(edg!=-1) {
      INTEGER(res)[edg + i * nedges]=1;
      edg=cnnctngEdg(edg, edge1, edge2, nedges);
    }
  }
  UNPROTECT(1);
  return res;
}
