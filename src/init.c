// generated using
// tools::package_native_routine_registration_skeleton()

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP ddirichlet_log_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP ddirichlet_log_vector(SEXP, SEXP, SEXP);
extern SEXP rdirichlet_matrix(SEXP, SEXP, SEXP);
extern SEXP rdirichlet_vector(SEXP, SEXP);
extern SEXP wght_LL_grad_alternative(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wght_LL_grad_common(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ddirichlet_log_matrix",    (DL_FUNC) &ddirichlet_log_matrix,     4},
    {"ddirichlet_log_vector",    (DL_FUNC) &ddirichlet_log_vector,     3},
    {"rdirichlet_matrix",        (DL_FUNC) &rdirichlet_matrix,         3},
    {"rdirichlet_vector",        (DL_FUNC) &rdirichlet_vector,         2},
    {"wght_LL_grad_alternative", (DL_FUNC) &wght_LL_grad_alternative, 15},
    {"wght_LL_grad_common",      (DL_FUNC) &wght_LL_grad_common,      10},
    {NULL, NULL, 0}
};

void R_init_DirichletReg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
