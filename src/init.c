#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_SUB(cvbmroc)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_SUB(cvbmrocpartial)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_SUB(cvbmrocfpf2tpf)(void *, void *, void *, void *, void *);
extern void F77_SUB(cvbmroctpf2fpf)(void *, void *, void *, void *, void *);
extern void F77_SUB(pbmroc)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_SUB(pbmrocpartial)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_SUB(pbmrocfpf2tpf)(void *, void *, void *, void *, void *);
extern void F77_SUB(pbmroctpf2fpf)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"cvbmroc", (DL_FUNC) &F77_SUB(cvbmroc), 8},
    {"cvbmrocpartial", (DL_FUNC) &F77_SUB(cvbmrocpartial), 7},
    {"cvbmrocfpf2tpf", (DL_FUNC) &F77_SUB(cvbmrocfpf2tpf), 5},
    {"cvbmroctpf2fpf", (DL_FUNC) &F77_SUB(cvbmroctpf2fpf), 5},
    {"pbmroc", (DL_FUNC) &F77_SUB(pbmroc), 8},
    {"pbmrocpartial", (DL_FUNC) &F77_SUB(pbmrocpartial), 7},
    {"pbmrocfpf2tpf", (DL_FUNC) &F77_SUB(pbmrocfpf2tpf), 5},
    {"pbmroctpf2fpf", (DL_FUNC) &F77_SUB(pbmroctpf2fpf), 5},
    {NULL, NULL, 0}
};

void R_init_MRMCaov(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
