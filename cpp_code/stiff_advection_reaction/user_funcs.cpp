#include "user_funcs.h"
#include "tensors.h"
#include "CONSTANTS.h"
#include "math.h"

void AuxInit( const dTensor2& xpts, dTensor2& aux )
{

    const int mpts = xpts.get_size(1);
    const int maux = aux.get_size(2);
    for( int i=1; i <= mpts; i++ )
    {

        const double x = xpts.get(i,1);
        aux.set(i, 1, sin(2.0 * PI*x ) );

    }


}

void Qinit( const dTensor2& xpts, dTensor2& qpts )
{

    const int mpts = xpts.get_size(1);
    const int meqn = qpts.get_size(2);

    const double width = 0.2;
    for( int i=1; i <= mpts; i++ )
    {

        const double x = xpts.get(i,1);
//      qpts.set(i, 1, sin(2.0 * PI*x ) );

        // cosine bump function.  this function has 5 derivatives //
        if( x > 0.5 && x < 0.7 )
        {
            double tmp = 1.0 + pow( cos( PI*(x-0.4) / width ), 6 );
            qpts.set(i, 1, tmp );
        }
        else
        {
            qpts.set(i,1,1.0);
        }

    }

}


void flux_func( double x, const dTensor1& q, dTensor1& f )
{

    // TODO - this should be reading the speed!
    f.set(1, 1.0 * q.get(1) );

}

void source_func( double x, const dTensor1& qpts, 
        const dTensor1& au, dTensor1& ps)
{

    //source term function
    const int meqn = qpts.get_size();

    double q1 = qpts.get(1);
    double f0 = au.get(1);

    // stiffness parameter
    double eps = 1e-2;
    ps.set(1, -1.0 / eps * ( q1 - f0 ) );

}

void flux_func_imp( double x, const dTensor1& qx, dTensor1& f )
{
    f.set(1, 0. );
}
