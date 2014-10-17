#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "tensors.h"
#include "CONSTANTS.h"

void func( const dTensor2& xpts, dTensor2& qpts )
{

    const int mpts = xpts.get_size(1);
    for( int i=1; i <= mpts; i++ )
    {

        const double x = xpts.get(i,1);
        qpts.set(i, 1, 1.0 + sin(10.0*PI*x ) );

    }

}

int main( int argc, char** argv )
{
    
    void integrate_on_cells( int sorder, 
        const dTensor2& node, dTensorBC2& q,  
        void (*Func)( const dTensor2& xpts, dTensor2& qpts ) );

    const int sorder = 2;

    int mx = 10;
    if( argc > 1 )
    {
        mx = atoi( argv[1] );
    }

    // grid spacing
    const double xlow   = 0.0;
    const double xhigh  = 1.0;
    const double dx = (xhigh-xlow) / mx;

    const int mbc  = 2;
    const int meqn = 1;

    // uniformly spaced node
    dTensor2 node( mx+1, 1 );
    for( int n=1; n <= mx+1; n++ )
    {
        node.set(n, 1, (n-1)*dx );
    }

    dTensorBC2 q(mx, meqn, mbc );
    integrate_on_cells( sorder, node, q, &func );

    for( int i=1; i <= mx; i++ )
    {
        printf("%2.8e\n", q.get(i,1) );
    }

    printf("  %d = q.get_size(1)\n", q.get_size(1) );
    return 0;
}
