#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _STAREG_em_lfdr(void *, void *, void *, void *);

// static const R_CallMethodDef CallEntries[] = {
//     {"_STAREG_em_lfdr", (DL_FUNC) &_STAREG_em_lfdr, 4},
//     {NULL, NULL, 0}
// };

// void R_init_STAREG(DllInfo *dll)
// {
//     R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
//     R_useDynamicSymbols(dll, FALSE);
// }

