#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_faltld.h"
#include "fxt_faltld_loc.h"

#include "fxt_math.h"
#include "fxt_fxtld.h"
#include "fxt_lagld.h"
#include "fxt_sarld.h"

#include "fxt_matll.h"
#include "fxt_vecll.h"
#include "fxt_matld.h"
#include "fxt_vecld.h"

#include "fxt_file.h"

/********************************************************************* 
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  generate preprocessed file for Fast Associated Legendre Transform
*********************************************************************/

static double res_prec;
static long res_mflop;

static void makeinfofile_f77(long p, long n, fxt_vecl *mv,
			     double prec, char *fname);

static long i_size;
static long d_size;
static long w_size;

static void create_legendre(int odd, fxt_vecll *x, fxt_vecll *w,
			    long n, long m, double prec, FILE *fout) {
  fxt_matll *lg;		/* the Legendre matrix */
  fxt_vecll *xx, *ww;		/* halved vectors */
  fxt_matld *dlg;		/* double-valued matrix */
  fxtld *fxt;			/* FXT structure */
  fxt_lagld *lag;		/* linear algorithm graph */
  fxt_sarld *ar;		/* simple array representation */
  fxt_vecld *dw;
  double err;

  /* save m */
  if (!odd)
    fxt_file_writelongs(fout, &m, 1);

  /* skip odd part for the 'last' order */
  if (odd && m == n)
    return;

  /* get Legendre function matrix */
  lg = fxt_legendre_matll(x, m, n, odd);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* get half of the Gaussian points */
  xx = fxt_vecll_clone(x, 0, lg->nrow - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* functions are polynomials of x**2 */
  fxt_vecll_sq(xx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* negate for upward sorting */
  fxt_vecll_neg(xx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* get the halved Gaussian weight vector */
  ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* double the value of duplicated points */
  fxt_vecll_scale(ww, 0, w->n / 2 - 1, 2.0L);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* square-roots, being multiplied twice */
  fxt_vecll_sqrt(ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* scale the Legendre matrix */
  fxt_matll_scalerow(lg, ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* create fxt: the main part */
  fxt = fxtld_from_matll(lg, xx, prec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* get the double-valued matrix */
  dlg = fxt_matll_to_matld(lg);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* convert to linear algorithm graph */
  lag = fxtld_get_lagld(fxt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* optimize as for the precision */
  fxt_matld_lagld_opt(lag, dlg, prec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* check error */
  err = fxt_matld_lagld_error(lag, dlg, POWERPREC, 0.0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (err > prec) {
#if DISALLOWHIGHERROR
    if (err > prec + 1e-16 * lg->nrow) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_faltld_preproc: precision unattained\n");
      return;
    } else
#endif
      fxt_error_set(FXT_ERROR_WARN,
		    "fxt_faltld_preproc: precision unattained\n");
  }

  res_prec = err;
  res_mflop = fxt_lagld_evaluate_mflop(lag);

  /* convert to simple array representation */
  ar = fxt_sarld_from_lagld(lag);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* descale sarld */
  dw = fxt_vecll_to_vecld(ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecld_einv(dw);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_sarld_scalerow(ar, dw);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecld_del(dw);

  /* save the results */
  fxt_sarld_save(fout, ar);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* count sizes */
  i_size += ar->nrow + ar->ncol + ar->ntmp + 1
    + ar->p[ar->ncol + ar->ntmp];
  d_size += ar->p[ar->ncol + ar->ntmp];
  
  if (w_size < ar->nrow + ar->ncol + ar->ntmp)
    w_size = ar->nrow + ar->ncol + ar->ntmp;

  /* deallocate everything */
  fxt_sarld_del(ar);
  fxt_lagld_del(lag);
  fxt_matld_del(dlg);
  fxtld_del(fxt);
  fxt_vecll_del(ww);
  fxt_vecll_del(xx);
  fxt_matll_del(lg);
}

void fxt_faltld_preproc(long p, long n, fxt_vecl *mv, double prec,
			char *fname) {
  long i;
  fxt_vecll *x, *w;
  fxt_vecld *v;
  long mflop0 = 0, dflop0 = 0;
  double errmax0 = 0.0;
  long data[3];
  FILE *fout;

  /* check null pointers */
  if (mv == NULL || fname == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_preproc: null pointer\n");
    return;
  }

  d_size = i_size = w_size = 0;

  /* check input */
  if (p < n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_genfaltld: irregal degree n\n");
    return;
  }

  /* open file */
  fout = fopen(fname, "wb");
  if (fout == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_genfaltld: cannot open %s\n", fname);
    return;
  }

  /* write header */
  fxt_file_writestr(fout, fxt_faltld_header);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* save size */
  data[0] = p;
  data[1] = n;
  data[2] = mv->n;

  fxt_file_writelongs(fout, data, 3);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* save precision */
  fxt_file_writedoubles(fout, &prec, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* create Gaussian points */
  fxt_gauss_vecll(p, &x, &w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* save Gaussian points */
  v = fxt_vecll_to_vecld(x);

  fxt_vecld_save(fout, v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecld_del(v);

  /* save Gaussian weights */
  v = fxt_vecll_to_vecld(w);

  fxt_vecld_save(fout, v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecld_del(v);

  for (i=0; i< mv->n; i++) {
    long m = mv->v[i];
    long mflop1 = 0, dflop1 = 0;
    double errmax1 = 0.0;

    /* check m */
    if (n < m) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_genflt: irregal order m\n");
      return;
    }

    printf("- Spherical(%ld %ld) %ld points:", n, m, p);
    fflush(stdout);

    /* create even half */
    create_legendre(0, x, w, n, m, prec, fout);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    mflop1 += res_mflop;
    if (errmax1 < res_prec)
      errmax1 = res_prec;

    /* create odd half */
    create_legendre(1, x, w, n, m, prec, fout);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    mflop1 += res_mflop;
    if (errmax1 < res_prec)
      errmax1 = res_prec;

    /* total results */
    dflop1 = (n - m + 1) * p / 2;

    printf(" error=%9.2e speedup=%4.2f\n", errmax1,
	   (double) dflop1 / (double) mflop1);

    mflop0 += mflop1;
    if (errmax0 < errmax1)
      errmax0 = errmax1;

    dflop0 += dflop1;
  }

  printf("- TOTAL max_error=%9.2e (req_prec=%9.2e) speedup=%4.2f\n",
	 errmax0, prec, (double) dflop0 / (double) mflop0);

  /* finalize file */
  fxt_file_writestr(fout, fxt_faltld_header);

  fxt_vecll_del(w);
  fxt_vecll_del(x);

  /* close file */
  fclose(fout);

  /* make data file for FORTRAN77 */
  makeinfofile_f77(p, n, mv, prec, fname);
}


static void makeinfofile_f77(long p, long n, fxt_vecl *mv,
			     double prec, char *fname) {
  FILE *f;  char *s;
  long mmax, nt, i;

  /* allocate info file name string */
  s = (char*) malloc(sizeof(char) * (strlen(fname) + 5));
  if (s == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_genfaltld: allocation failed\n");
    return;
  }

  /* make info file name */
  strcpy(s, fname);
  strcat(s, "_f77");

  /* open info file */
  f = fopen(s, "w");
  if (f == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_genfaltld: cannot open %s\n", s);
    return;
  }

  free(s);

  /* determine max m and number of transforms */
  mmax = nt = 0;
  for (i=0; i< mv->n; i++) {
    if (mmax < mv->v[i])
      mmax = mv->v[i];

    if (mv->v[i] == n)
      nt += 1;
    else
      nt += 2;
  }

  fprintf(f, "C ========= FXTPACK information file for FORTRAN77\n");
  fprintf(f, "C --------- data file name\n");
  fprintf(f, "      character(*) fxtfnm\n");
  fprintf(f, "      parameter (fxtfnm = '%s')\n", fname);
  fprintf(f, "\n");
  fprintf(f, "C --------- version string\n");
  fprintf(f, "      character(*) fxtfhd\n");
  fprintf(f, "      parameter (fxtfhd = '%s')\n", fxt_faltld_header);
  fprintf(f, "\n");
  fprintf(f, "C --------- length of version string\n");
  fprintf(f, "      integer fxthdl\n");
  fprintf(f, "      parameter (fxthdl = %ld)\n",
	  (long) strlen(fxt_faltld_header));
  fprintf(f, "\n");
  fprintf(f, "C --------- number of points\n");
  fprintf(f, "      integer fxtp\n");
  fprintf(f, "      parameter (fxtp = %ld)\n", p);
  fprintf(f, "\n");
  fprintf(f, "C --------- maximum degree n\n");
  fprintf(f, "      integer fxtnmx\n");
  fprintf(f, "      parameter (fxtnmx = %ld)\n", n);
  fprintf(f, "\n");
  fprintf(f, "C --------- maximum order  m\n");
  fprintf(f, "      integer fxtmmx\n");
  fprintf(f, "      parameter (fxtmmx = %ld)\n", mmax);
  fprintf(f, "\n");
  fprintf(f, "C --------- number of half-transforms\n");
  fprintf(f, "      integer fxtnt\n");
  fprintf(f, "      parameter (fxtnt = %ld)\n", nt);
  fprintf(f, "\n");
  fprintf(f, "C --------- integer data size\n");
  fprintf(f, "      integer fxtisz\n");
  fprintf(f, "      parameter (fxtisz = %ld)\n", i_size);
  fprintf(f, "\n");
  fprintf(f, "C --------- double precision data size\n");
  fprintf(f, "      integer fxtdsz\n");
  fprintf(f, "      parameter (fxtdsz = %ld)\n", d_size);
  fprintf(f, "\n");
  fprintf(f, "C --------- precision of the fast transform\n");
  fprintf(f, "      double precision fxtprc\n");
  fprintf(f, "      parameter (fxtprc = %e)\n", prec);
  fprintf(f, "\n");
  fprintf(f, "C --------- maximum working array size\n");
  fprintf(f, "      integer fxtwmx\n");
  fprintf(f, "      parameter (fxtwmx = %ld)\n", w_size);

  fclose(f);
}
