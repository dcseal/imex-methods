#ifndef _MATRIX_SOLVE_H_
#define _MATRIX_SOLVE_H_

#include <stdio.h>
#include <stdlib.h>
#include "tensors.h"
#include "petscksp.h"

typedef struct
{

    Mat A;    /* linear system matrix */
    KSP ksp;  /* linear solver context */
    PC  pc;   /* preconditioner context */

    Vec x;    /* vector of unknowns */
    Vec b;    /* right hand side vector */

} matrix_objs;

int matrix_solve( matrix_objs& objs, double eps, const int mbc_order,
    double dx, const dTensorBC2& f, dTensorBC2& q );

int InitObjs( matrix_objs &objs, int mx, double dx, double as );
int CleanObjs( matrix_objs & objs );

#endif
