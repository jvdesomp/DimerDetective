/*
 * dd_thal.c - thin wrapper around primer3's thal() for DimerDetective.
 *
 * Reproduces primer3-py's ThermoAnalysis.calc_end_stability() (alignment
 * type thal_end1, mode THL_GENERAL, dimer = 1) and additionally asks thal()
 * for the secondary structure string. With print_output = 1, thal() calls
 * drawDimer() instead of calcDimer(); both use the same formula for dG, so
 * the returned dG is identical to calc_end_stability() and the structure is
 * the same SEQ/STR block that `ntthal -a END1` prints.
 *
 * The thermodynamic parameters are read from the primer3_config directory,
 * exactly as primer3-py does at import time.
 *
 * Copyright (C) 2026 Ghent University. Licensed under GPL-2.0, like primer3.
 */

#include <stdlib.h>
#include <string.h>

#include "primer3/thal.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define DD_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define DD_EXPORT
#endif

static thal_results last;
static char *struct_buf = NULL;
static int loaded = 0;

/* Load parameters from `path` (must end with a slash). Returns 0 on success. */
DD_EXPORT int dd_init(const char *path)
{
  thal_parameters tp;
  thal_results o;
  int rc;

  if (loaded) return 0;
  memset(&o, 0, sizeof(o));
  thal_set_null_parameters(&tp);
  if (thal_load_parameters(path, &tp, &o) != 0) return -1;
  rc = get_thermodynamic_values(&tp, &o);
  thal_free_parameters(&tp);
  if (rc != 0) return -2;
  loaded = 1;
  return 0;
}

/*
 * 3' end stability of seq1 against seq2 (both 5'->3', upper case ACGT).
 * Concentrations as in primer3-py: mv, dv, dntp in mM, dna_conc in nM,
 * temp_c in degrees Celsius. Returns dG in cal/mol (0 if no structure).
 * The structure and the other values are available through dd_last_*().
 */
DD_EXPORT double dd_end_stability(const char *seq1, const char *seq2,
                                  double mv, double dv, double dntp,
                                  double dna_conc, double temp_c,
                                  int max_loop)
{
  thal_args a;
  size_t n = strlen(seq1) + strlen(seq2);

  memset(&a, 0, sizeof(a));
  a.type = thal_end1;
  a.maxLoop = max_loop;
  a.mv = mv;
  a.dv = dv;
  a.dntp = dntp;
  a.dna_conc = dna_conc;
  a.temp = temp_c + 273.15;
  a.dimer = 1;

  free(struct_buf);
  struct_buf = (char *) malloc(n * 4 + 64);
  if (struct_buf == NULL) return 0.0;
  struct_buf[0] = '\0';

  memset(&last, 0, sizeof(last));
  last.sec_struct = struct_buf;
  thal((const unsigned char *) seq1, (const unsigned char *) seq2,
       &a, THL_GENERAL, &last, 1);
  last.sec_struct = NULL;
  if (last.msg[0] != '\0' || last.no_structure) {
    struct_buf[0] = '\0';
  }
  return last.dg;
}

DD_EXPORT const char *dd_last_structure(void) { return struct_buf ? struct_buf : ""; }
DD_EXPORT const char *dd_last_message(void) { return last.msg; }
DD_EXPORT int dd_last_no_structure(void) { return last.no_structure; }
DD_EXPORT double dd_last_tm(void) { return last.temp; }
DD_EXPORT double dd_last_dh(void) { return last.dh; }
DD_EXPORT double dd_last_ds(void) { return last.ds; }
