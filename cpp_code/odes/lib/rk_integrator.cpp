#include <stdio.h>
#include <assert.h>
#include "defs.h"
#include "req_state.h"
#include "rk_coeffs.h"

void construct_l( const dTensor2* node, dTensorBC2* q, dTensorBC2* aux, dTensorBC2* L );

// This function steps q from qn to qnp1.
// This function doesn't attempt to save q by any means.
// The choice of integration technique is given by the rk_coeffs.h file.
//
// This routine could be much more efficient if taylored to work with the
// butcher tableau rather than brute forcing its way through the tableau.
//
void rk_integrator( req_state* top_state, 
                        dTensorBC2** kE, dTensorBC2& yi )
{

    dTensorBC2* q   = top_state->get_q();
    dTensorBC2* aux = top_state->get_aux();

    const double dt = top_state->get_dt();

    const int mx   = q->get_size(1);
    const int meqn = q->get_size(2);
    const int mbc  = q->get_mbc();


    // Run the s-stage method //
    for( int i=0; i < s; i++ )
    {

        #pragma omp parallel for
        for( int ix=1; ix <= mx; ix++ )
        for( int m=1; m <= meqn; m++ )
        {

            double tmp = q->get(ix,m);
            for( int j=0; j < i; j++ )
            {
                tmp += ( dt * AE[i][j] ) * kE[j]->get(ix,m);

                // TODO implicit part too!
                // tmp += dt * AI[i,j] * kI[j]->get(ix,m)

            }
            yi.set( ix, m, tmp );

        }

        // TODO -- implicit solve needs to happen here!!
        // yi = yi + dt A[i,i] * fI( yi );

        // evaluate the right hand side:
        construct_l( top_state->get_node(), &yi, aux, kE[i] );

        // This should really look like:
        // yi->construct_l( kE[i] );

    }

    // update the solution:
    #pragma omp parallel for
    for( int ix=1; ix <= mx; ix++ )
    for( int m=1; m <= meqn; m++ )
    {

        double tmp = q->get(ix,m);
        for( int n=0; n < s; n++ )
        {
            tmp+= dt * bE[n] * kE[n]->get( ix, m );
        }

        q->set( ix, m, tmp );

    }


/*
 * Classical RK4 Integrator
 *

    // stage 1:
    construct_l( xlow, dx, q, kE[0] );

    // stage 2:
    for( int ix=1; ix <= mx; ix++ )
    for( int m=1; m <= meqn; m++ )
    {
        double tmp = q->get(ix,m) + 0.5 * dt * kE[0]->get(ix,m);
        yi.set(ix,m, tmp );
        construct_l( xlow, dx, &yi, kE[1] );
    }

    // stage 3:
    for( int ix=1; ix <= mx; ix++ )
    for( int m=1; m <= meqn; m++ )
    {
        double tmp = q->get(ix,m) + 0.5 * dt * kE[1]->get(ix,m);
        yi.set(ix,m, tmp );
        construct_l( xlow, dx, &yi, kE[2] );
    }

    // stage 4:
    for( int ix=1; ix <= mx; ix++ )
    for( int m=1; m <= meqn; m++ )
    {
        double tmp = q->get(ix,m) + dt * kE[2]->get(ix,m);
        yi.set(ix,m, tmp );
        construct_l( xlow, dx, &yi, kE[3] );
    }

*/

}
