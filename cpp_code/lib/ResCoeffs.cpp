#include "tensors.h"
#include <stdlib.h>
#include <stdio.h>

/*
 * Function for computing the coefficients necessary for integrating the
 * residual.  Integration on the 'i'-th sub-interval is of the form:
 *
 *    N_i = coeffs(i,:) * Fvals(:);
 *
 * where Fvals(:) are the function values at each of the num_points.
 *
 */
void ResCoeffs( dTensor2& coeffs )
{

    // For the most efficient (forward/backward euler), 
    //num_points = num_intervals + 1.
    const int num_intervals = coeffs.get_size(1);
    const int num_points    = coeffs.get_size(2);


    switch( num_points )
    {

        case 2:
    
            coeffs.set(1, 1, 1.0);
            coeffs.set(1, 2, 1.0);
            break;

        case 3:

            coeffs.set(1,1,  5.0 / 12.0 );
            coeffs.set(1,2,  2.0 /  3.0 );
            coeffs.set(1,3, -1.0 / 12.0 );

            coeffs.set(2,1, -1.0 / 12.0 );
            coeffs.set(2,2,  2.0 /  3.0 );
            coeffs.set(2,3,  5.0 / 12.0 );
            break;

        case 4:

            coeffs.set(1,1, 0.25 );
            coeffs.set(1,2, 19.0 / 36.0 );
            coeffs.set(1,3, -5.0 / 36.0  );
            coeffs.set(1,4,  1.0 / 36.0 );
            
            coeffs.set(2,1, -coeffs.get(1,4) );
            coeffs.set(2,2, 13.0 / 36.0     );
            coeffs.set(2,3, coeffs.get(2,2) );
            coeffs.set(2,4, coeffs.get(2,1) );

            coeffs.set(3,1, coeffs.get(1,4) );
            coeffs.set(3,2, coeffs.get(1,3) );
            coeffs.set(3,3, coeffs.get(1,2) );
            coeffs.set(3,4, coeffs.get(1,1) );
            break;

    }

}

/*
 * Function for integrating the residual.  Res should be function or Q values at
 * each of the intermediate grid points.  N is integrated values of Res for each
 * of the sub-intervals.
 *
 * TODO - It would be nice to formulate this entire thing in terms of a MOL
 * discretization, rather than expect the users data to conform to a specific
 * BC[n] type.
 *
 */
void IntegrateResidual( const dTensorBC3& Res, dTensorBC3& N)
{

    // space-time grid data //
    const int mx            = N.get_size(1);
    const int meqn          = N.get_size(2);
    const int num_intervals = N.get_size(3);
    const int mbc           = N.get_mbc();

    const int num_points = num_intervals+1;

    // Grab the necessary coefficients
    dTensor2 coeffs( num_points-1, num_points );
    ResCoeffs( coeffs );
    
    for( int i=1-mbc; i <= mx + mbc; i++ )
    for( int m=1; m <= meqn; m++ )
    for( int n=1; n <= num_intervals; n++ )
    {
        double tmp = 0.0;
        for( int k=1; k <= num_points; k++ )
        {
            tmp += coeffs.get(n,k) * Res.get(i,m,k);
        }
        N.set(i,m,n, tmp );

    }

}
