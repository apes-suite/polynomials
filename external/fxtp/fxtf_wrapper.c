// Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
// Copyright (c) 2015 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
//
// Parts of this file were written by Harald Klimach and Kay Langhammer
// for University of Siegen.
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

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
