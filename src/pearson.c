#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "pearson.h"

static const R_CMethodDef CEntries[] = {
    {"FPUcheck", (DL_FUNC) &FPUcheck, 0},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"logPearsonIVnorm", (DL_FUNC) &logPearsonIVnorm, 3},
    {"PearsonIVnorm",    (DL_FUNC) &PearsonIVnorm,    3},
    {"rPearsonIVk",      (DL_FUNC) &rPearsonIVk,      6},
    {"rPearsonIVlogK",   (DL_FUNC) &rPearsonIVlogK,   6},
    {"F21D",             (DL_FUNC) &F21D,             6},
    {"F21Da1bR",         (DL_FUNC) &F21Da1bR,         6},
    {"F21Da1cR",         (DL_FUNC) &F21Da1cR,         6},
    {"F21DaR",           (DL_FUNC) &F21DaR,           6},
    {"F21DD",            (DL_FUNC) &F21DD,            6},
    {"F21DDa1bR",        (DL_FUNC) &F21DDa1bR,        6},
    {"F21DDa1cR",        (DL_FUNC) &F21DDa1cR,        6},
    {"F21DDaR",          (DL_FUNC) &F21DDaR,          6},
    {"F21QD",            (DL_FUNC) &F21QD,            6},
    {"F21QDa1bR",        (DL_FUNC) &F21QDa1bR,        6},
    {"F21QDa1cR",        (DL_FUNC) &F21QDa1cR,        6},
    {"F21QDaR",          (DL_FUNC) &F21QDaR,          6},
    {NULL, NULL, 0}
};

void R_init_PearsonDS(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
