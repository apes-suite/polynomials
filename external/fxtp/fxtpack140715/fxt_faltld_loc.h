#ifndef __FXT_FALTLD_LOC_H_INCLUDED__
#define __FXT_FALTLD_LOC_H_INCLUDED__

#include "fxt_faltld.h"
#include "fxt_vecld.h"
#include "fxt_sarld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  internal header file for Fast Associated Legendre Transform
*********************************************************************/

/*** header string for the file ***/
extern const char* fxt_faltld_header;

/*** data structure for each associated Legendre transform ***/
struct FALTLD_ALT {
  fxt_faltld *par;		/* parent structure */

  long m;			/* order */
  
  fxt_sarld *ear;		/* sarld for even part */
  fxt_sarld *oar;		/* sarld for odd  part */
};

#endif
