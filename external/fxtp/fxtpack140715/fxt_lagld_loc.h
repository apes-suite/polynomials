#ifndef __FXT_LAGLD_LOC_H__
#define __FXT_LAGLD_LOC_H__

#include "fxt_lagld.h"

/*********************************************************************
   inner header file for linear algorithm graph
*********************************************************************/

typedef struct FXT_LAGLD_ENT fxt_lagld_ent;

/*** the LAGLD data structure ***/
struct FXT_LAGLD {
  long ncol;			/* size of the input vector */
  long nrow;			/* size of the output vector */

  fxt_lagld_vec *ivec;		/* entries for the input vector */
  fxt_lagld_vec *ovec;		/* entries for the output vector */

  fxt_lagld_var *valloc;	/* allocated VAR's */
  fxt_lagld_ent *ealloc;	/* allocated ENT's */

  long na_var, nu_var;		/* allocated/used VAR's */
  long na_ent, nu_ent;		/* allocated/used ENT's */
};


/*** variable data structure ***/
struct FXT_LAGLD_VAR {
  int io;			/* 1 for input, 2 for output */

  long n_in;			/* number of inputs */
  long n_out;			/* number of outputs */

  fxt_lagld_ent *ient;		/* list of input entries */
  fxt_lagld_ent *oent;		/* list of output entries */

  double v;			/* the value */
  int c;			/* counter */

  double iw;			/* weight from inputs */
  double ow;			/* weight from outputs */
};


/*** matrix element ***/
struct FXT_LAGLD_ENT {
  double v;			/* the value */
  double w;			/* weight */

  fxt_lagld_var *ivar;		/* the input variable */
  fxt_lagld_var *ovar;		/* the output variable */

  fxt_lagld_ent *inext;		/* next for ient list */
  fxt_lagld_ent *onext;		/* next for oent list */
};

void fxt_lagld_setweight(fxt_lagld *lag);

#endif
