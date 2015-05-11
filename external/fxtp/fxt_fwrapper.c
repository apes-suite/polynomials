#include "fxt_vecld.h"
#include "fxt_flptld.h"

void fxtf_flptld_evl(double *v, int vlen, fxt_flptld *flpt,
                     double *u, int ulen, fxt_vecld *w) {
  fxt_vecld vvec;
  fxt_vecld uvec;

  vvec.n = vlen;
  vvec.v = v;
  uvec.n = ulen;
  uvec.v = u;

  fxt_flptld_evl(&vvec, flpt, &uvec, w);
}

void fxtf_flptld_exp(double *u, int vlen, fxt_flptld *flpt,
                     double *v, int ulen, fxt_vecld *w) {
  fxt_vecld vvec;
  fxt_vecld uvec;

  vvec.n = vlen;
  vvec.v = v;
  uvec.n = ulen;
  uvec.v = u;

  fxt_flptld_exp(&uvec, flpt, &vvec, w);
}
