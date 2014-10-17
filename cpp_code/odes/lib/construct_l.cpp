#include <stdio.h>
#include "defs.h"

#include "req_state.h"
#include "user_funcs.h"
void set_boundary( dTensorBC2& q );

void construct_l( const dTensor2* node, dTensorBC2* q, 
    dTensorBC2* aux, dTensorBC2* L )
{

    const int mx   = q->get_size(1);
    const int meqn = q->get_size(2);
    const int maux = aux->get_size(2);
    const int mbc  = q->get_mbc();
    set_boundary( *q );


    const double xlow = node->get(1,1);
    const double dx   = node->get(2,1) - node->get(1,1);

    dTensorBC2* F   = new dTensorBC2(mx+1, meqn, mbc );

    // stencil
    dTensorBC1* stencil = new dTensorBC1( 1, mbc );
    switch( mbc )
    {

        case 1:

            // centered stencil (again!)
            stencil->set(0, 0.5);
            stencil->set(1, 0.5);

        case 2:

            // centered stencil
            stencil->set( -1, -1.0/12.0 );
            stencil->set(  0,  7.0/12.0 );
            stencil->set(  1,  7.0/12.0 );
            stencil->set(  2, -1.0/12.0 );
    }

    #pragma omp parallel for
    for( int i=1; i <= mx+1; i++ )
    { 

        // sample the function at the cell interface
        dTensor1 qimh( meqn );
        for( int m=1; m <= meqn; m++ )
        {
            double qi = 0.0;
            for( int l=-mbc; l < mbc; l++ )
            {
                qi += stencil->get(l+1) * q->get( i+l, m );
            }

            qimh.set(m, qi );
        }

        // evaluate the flux function at the interface
        dTensor1 fval( meqn );
        flux_func( node->get(i,1), qimh, fval );

        // save the result
        for( int m=1; m <= meqn; m++ )
        { F->set( i, m, fval.get(m) ); }

    }

    #pragma omp parallel for
    for( int i=1; i <= mx; i++ )
    {

        // Grab the source term
        dTensor1 ps( meqn );
        dTensor1 qp( meqn );
        dTensor1 au( maux );
        for( int m=1; m <= meqn; m++ )
        {
            ps.set(m, 0. );
            qp.set(m, q->get(i,m) );
        }
        for( int m=1; m <= maux; m++ )
        {
            au.set(m, aux->get(i,m) );
        }
        source_func( 0.5*(node->get(i+1,1) - node->get(i,1) ), qp, au, ps );
    
        for( int m=1; m <= meqn; m++ )
        {
            // flux term
            double fl = -( F->get(i+1,m) - F->get(i,m) ) / dx;

            L->set( i, m, fl + ps.get(m) );
        }

    }

    delete F;
    delete stencil;

}
