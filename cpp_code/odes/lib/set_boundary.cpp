#include "tensors.h"

// periodic boundary conditions here:
void set_boundary( dTensorBC2& q )
{

    const int mx   = q.get_size(1);
    const int meqn = q.get_size(2);
    const int mbc  = q.get_mbc();

    // TODO!  INCORPORATE BOUNDARY CONDITIONS OTHER THAN PERIODIC!
    int mbc_type [] = {0, 0};

    if( mbc_type[0] == 0 )
    { 
        for( int i=1; i <= mbc; i++ )
        for( int m=1; m <= meqn; m++ )
        {
            // right boundary cells:
            q.set( mx+i, m, q.get( i, m ) );

            // left boundary cells:
            q.set( 1-i,   m, q.get( (mx+1)-i, m ) );

        }
        return;
    }

    // dirichlet on left hand side:
    if( mbc_type[0] == 1 )
    {

        // left hand boundary values
        dTensor1 A(meqn);
        for( int m=1; m <= meqn; m++ )
        { A.set(m,0.0); }

        for( int m=1; m <= meqn; m++ )
        {


            // left hand boundary values (TODO - CALL A FUNCTION TO SET THIS?)
            if( mbc == 1 )
            {
                q.set(0, m, 2.0*A.get(m) - q.get(1,m) );
            } 
            else if( mbc == 2 )
            {
                q.set(0, m, 
                    4.0*A.get(m) - 13.0/3.0*q.get(1,m) +
                    5.0/3.0*q.get(2,m)-1.0/3.0*q.get(3,m) );
                q.set(-1, m, 
                    -12.0*A.get(m) + 7.0*( q.get(1,m) + q.get(0,m) )
                    - q.get(2,m) );
            }

        }
    }

    // dirichlet on right hand side:
    if( mbc_type[1] == 1 )
    {

        // right hand boundary values
        dTensor1 B(meqn);
        for( int m=1; m <= meqn; m++ )
        { B.set(m,0.0); } 

        for( int m=1; m <= meqn; m++ )
        {


            if( mbc == 1 )
            {
                q.set(mx+1, m, 2.0*B.get(m) - q.get(mx,1) );
            } 
            else if( mbc == 2 )
            {
                q.set(mx+1, m, 
                    4.0*B.get(m) + 1.0/3.0*( -13.0*q.get(mx,m) +
                    5.0*q.get(mx-1,m) - q.get(mx-2,m) ) );

                q.set(mx+2, m, 
                    4.0*B.get(m) + 1.0/3.0*( -13.0*q.get(mx+1,m) +
                    5.0*q.get(mx,m) - q.get(mx-1,m) ) );

            }

        }
    }

    // outflow on right hand side
    if( mbc_type[1] == 2 )
    {

        for( int m=1; m <= meqn; m++ )
        {

            if( mbc == 1 )
            {
                q.set(mx+1, m, 2.0*q.get(mx,1) - q.get(mx-1,m) );
            } 
            else if( mbc == 2 )
            {

                double tmp1 = -q.get(mx-3,m)+4.0*q.get(mx-2,m)
                              -6.0*q.get(mx-1,m)+4.0*q.get(mx,m);

                double tmp2 = -q.get(mx-2,m)+4.0*q.get(mx-1,m)
                              -6.0*q.get(mx,m)+4.0*q.get(1,m);

                q.set( mx+1, m, tmp1 );
                q.set( mx+2, m, tmp2 );

            }

        }
    }


}
