#include <stdlib.h>
#include <stdio.h>
#include "CONSTANTS.h"
#include "tensors.h"

// This function should be called exactly once
void set_gauss_points( dTensor1* w1d, dTensor1* x1d, int mpoints )
{

    switch( mpoints )
    {
        
        case 1:
            w1d->set(1, 2.0 );
            x1d->set(1, 0.0 );
            
            break;

        case 2:
            x1d->set(1, -1/sq3);
            x1d->set(2, 1/sq3);

            w1d->set(1, 1.0 );
            w1d->set(2, 1.0 );

            break;


        case 3:

            x1d->set(1, -sq3/sq5 );
            x1d->set(2,      0.0 );
            x1d->set(3,  sq3/sq5 );

            w1d->set(1, 5.0 / 9.0 );
            w1d->set(2, 8.0 / 9.0 );
            w1d->set(3, 5.0 / 9.0 );

            break;

        case 4:
            x1d->set(1, sqrt( ( 3.0 - 2.0*sq2*sq3/sq5 ) / 7.0 ) );
            x1d->set(2, - x1d->get(1) );

            x1d->set(3, sqrt( ( 3.0 + 2.0*sq2*sq3/sq5 ) / 7.0 ) );
            x1d->set(4, - x1d->get(3) );

            w1d->set(1, ( 18.0+sqrt( 30.0 ) ) / 36.0 );
            w1d->set(2, ( 18.0+sqrt( 30.0 ) ) / 36.0 );

            w1d->set(3, ( 18.0-sqrt( 30.0 ) ) / 36.0 );
            w1d->set(4, ( 18.0-sqrt( 30.0 ) ) / 36.0 );

            break;

        default:

            perror("Need to select a valid number of points\n");

            printf("  mpoitns = %d\n", mpoints );
    }

}

void integrate_on_cells( int sorder, 
        const dTensor2& node, dTensorBC2& q,  
        void (*Func)( const dTensor2& xpts, dTensor2& qpts ) 
    )
{

    const int mx   = q.get_size(1);
    const int meqn = q.get_size(2);
//  const int mbc  = q.get_mbc();

    // ?? TODO complete overkill here ... (TODO?)
    const int mpoints = 4;

    // weight and points
    dTensor1* x1d = new dTensor1( mpoints ); 
    dTensor1* w1d = new dTensor1( mpoints ); 
    set_gauss_points( w1d, x1d, mpoints );

    for( int i=1; i <= mx; i++ )
    {

        // values for q at each point

        dTensor2 qpts( x1d->get_size(), meqn );

        // cell center and quadrature points:
        const double b  = node.get(i+1,1);
        const double a  = node.get(i,1);
        const double xc = 0.5*( a + b );
        dTensor2 xpts( mpoints, 1 );
        for( int n=1; n <= mpoints; n++ )
        {
            xpts.set( n, 1, xc + 0.5*( (b-a)*x1d->get(n) ) );
        }

        // move the standard points to the new cell
        Func( xpts, qpts );
        for( int m=1; m <= meqn; m++ )
        {
            double qc = 0.0;
            for( int n=1; n <= mpoints; n++ )
            {
                qc += 0.5 * w1d->get(n) * qpts.get(n, m);
            }
            q.set(i, m, qc );
        }
    }

    delete x1d;
    delete w1d;

}
