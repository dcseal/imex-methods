#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "tensors.h"
#include "CONSTANTS.h"


void integrate_on_cells( int sorder, 
        const dTensor3& node, dTensorBC3& q,  
        void (*Func)( const dTensor2& xpts, dTensor2& qpts ) 
    );

void func( const dTensor2& xpts, dTensor2& qpts )
{

    const int mpts = xpts.get_size(1);
    for( int i=1; i <= mpts; i++ )
    {

        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        qpts.set(i, 1, sin(2.0*PI*x)*sin(2.0*PI*y) );

    }

}

int main( int argc, char** argv )
{

    for( int i=1; i<=5; i++ )
    {
        dTensorBC1 stencil(1,1);
        stencil.set(0,1.0);
        stencil.set(1,1.0);

    }

    const int mx   = 20;
    const int my   = 20;
    const int meqn = 1;
    const int mbc  = 2;


    const double xlow = 0.;  const double xhigh = 1.0; 
    const double ylow = 0.;  const double yhigh = 1.0; 

    const double dx = (xhigh-xlow) / mx;
    const double dy = (yhigh-ylow) / my;
    
    dTensor3 node(mx+1,my+1,2);
    for( int i=1; i <= mx+1; i++ )
    for( int j=1; j <= mx+1; j++ )
    {
        node.set(i,j, 1, xlow + (i-1)*dx );
        node.set(i,j, 2, ylow + (j-1)*dy );
    }

    dTensorBC3 q(mx,my,meqn,mbc);
    integrate_on_cells( 4, node, q, &func );

    dTensorBC3& qcp = q;

    printf("[ ");
    for( int i=1; i <= mx; i++ )
    for( int j=1; j <= my; j++ )
    for( int m=1; m <= meqn; m++ )
    {
        //printf("%2.3f %2.3f %2.5e\n", node.get(i,j,1), node.get(i,j,2),  q.get(i,j,m) );
        printf("%2.15e, ", qcp.get(i,j,m) );
    }
    printf("];\n");

  
    return 0;
}
