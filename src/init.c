#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void LC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PMD(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void LQuad(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void LifeTable(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void adjust_mx(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"LC", (DL_FUNC) &LC, 17},
    {"PMD", (DL_FUNC) &PMD, 19},
    {"LQuad", (DL_FUNC) &LQuad, 16},
    {"LifeTable", (DL_FUNC) &LifeTable, 12},
    {"adjust_mx", (DL_FUNC) &adjust_mx, 11},
    {NULL, NULL, 0}
};

void R_init_MortCast(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
