#include "user_funcs.h"
#include "tensors.h"
#include "CONSTANTS.h"
#include "math.h"

void AuxInit( const dTensor2& xpts, dTensor2& aux )
{

    const int mpts = xpts.get_size(1);
    const int maux = aux.get_size(2);

}

void Qinit( const dTensor2& xpts, dTensor2& qpts )
{

    const int mpts = xpts.get_size(1);
    const int meqn = qpts.get_size(2);

    const double width = 0.2;
    for( int i=1; i <= mpts; i++ )
    {

        const double x = xpts.get(i,1);
        qpts.set(i, 1, sin(2.0 * PI*x ) );

        // cosine bump function.  this function has 5 derivatives //
//      if( x > 0.5 && x < 0.7 )
//      {
//          double tmp = 1.0 + pow( cos( PI*(x-0.4) / width ), 6 );
//          qpts.set(i, 1, tmp );
//      }
//      else
//      {
//          qpts.set(i,1,1.0);
//      }



    }

}

void flux_func( double x, const dTensor1& q, dTensor1& f )
{

    // TODO - this should be reading the speed from a parameters file!
    // const double u = 1.0;
    const double u = 0.0;
    f.set(1, u * q.get(1) );

}

void flux_func_imp( double x, const dTensor1& qx, dTensor1& f )
{

    // TODO - this should be reading the speed from a parameters file!
    //const double eps = 1.0;
    const double eps = 1e-2;
    f.set(1, eps * qx.get(1) );

}

void source_func( double x, const dTensor1& q, const dTensor1& au, dTensor1& s )
{

    // TODO - this should be reading the speed from a parameters file!
    s.set(1, 0.0 );

}
