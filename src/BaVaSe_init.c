#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void eZSBF(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
	{"eZSBF",  (DL_FUNC) &eZSBF, 6},
		{NULL, NULL, 0}
};

void R_init_BaVaSe(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
