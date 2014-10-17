#ifndef _USER_FUNCS_
#define _USER_FUNCS_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tensors.h"

void Qinit( const dTensor2&  xpts, dTensor2& q );
void AuxInit( const dTensor2&  xpts, dTensor2& aux );
void FluxFunc( const dTensor2& xpts, dTensor2& q);
void flux_func( double x, const dTensor1& q, dTensor1& f );
void source_func( double x, const dTensor1& q, const dTensor1& au, dTensor1& f );

#endif
