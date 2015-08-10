/* Interface wrappers to allow interaction with Fortran.
   These fxtf routines are needed to encapsulate the fxt_vecld data structure
   and enable the direct passing of Fortran arrays to the FXTPACK routines.
 */
#include "fxt_vecld.h"
#include "fxt_flptld.h"
#include "fxt_faltld.h"

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

void fxtf_flptld_exp(double *u, int ulen, fxt_flptld *flpt,
                     double *v, int vlen, fxt_vecld *w) {
  fxt_vecld vvec;
  fxt_vecld uvec;

  vvec.n = vlen;
  vvec.v = v;
  uvec.n = ulen;
  uvec.v = u;

  fxt_flptld_exp(&uvec, flpt, &vvec, w);
}

void fxtf_faltld_evl(double *v, int vlen, fxt_faltld *falt, long m, 
                     double *u, int ulen, fxt_vecld *w) {
  fxt_vecld vvec;
  fxt_vecld uvec;

  vvec.n = vlen;
  vvec.v = v;
  uvec.n = ulen;
  uvec.v = u;

  fxt_faltld_evl(&vvec, falt, m, &uvec, w);
}


void fxtf_faltld_exp(double *u, int ulen, fxt_faltld *falt, long m, 
                     double *v, int vlen, fxt_vecld *w) {
  fxt_vecld uvec;
  fxt_vecld vvec;

  uvec.n = ulen;
  uvec.v = u;
  vvec.n = vlen;
  vvec.v = v;

  fxt_faltld_exp(&uvec, falt, m, &vvec, w);
}
