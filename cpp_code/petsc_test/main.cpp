#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "petscksp.h"
#include "tensors.h"
#include "CONSTANTS.h"

int matrix_solve( const int mbc_order, const dTensorBC2& qex, double dx, const dTensorBC2& f,
dTensorBC2& q, int argc, char** args);

void integrate_on_cells( int sorder, 
        const dTensor2& node, dTensorBC2& q,  
        void (*Func)( const dTensor2& xpts, dTensor2& qpts ) 
    );


void frhs( const dTensor2& xpts, dTensor2& fpts )
{
    const int mpts = fpts.get_size(1);
    const int meqn = fpts.get_size(2);
    for( int i=1; i <= mpts; i++ )
    {
        const double x = xpts.get(i,1);
        fpts.set(i, 1, -pow(2.0*PI, 2 ) * sin( 2.0 * PI * x ) );
    }

}

void Qexact( const dTensor2& xpts, dTensor2& Qpts )
{
    const int mpts = Qpts.get_size(1);
    const int meqn = Qpts.get_size(2);
    for( int i=1; i <= mpts; i++ )
    {
        const double x = xpts.get(i,1);
        Qpts.set(i, 1, -sin( 2.0 * PI * x ) );
    }

}

int main(int argc,char **args)
{

    const int mx        = 1000;
    const int mbc_order = 1;

    dTensorBC2 fbc(mx,1,2);
    dTensorBC2 qbc(mx,1,2);
    dTensorBC2 qex(mx,1,2);

    // set up grid
    const double xhigh = 1.0;  
    const double xlow  = 0.0;  
    const double dx    = (xhigh-xlow)/mx;

    dTensor2 node(mx+1,1);
    for( int i=1; i <= mx+1; i++ )
    { node.set(i, 1, xlow + (i-1)*dx ); } 
    integrate_on_cells( 4, node, fbc,  frhs   );
    integrate_on_cells( 4, node, qex,  Qexact );
    qbc.CopyFrom( &qex );


    char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
    PetscInitialize(&argc,&args,(char *)0,help);
    int result = matrix_solve(mbc_order, qex, dx, fbc, qbc, argc, args );
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);

    double norm = 0.0;
    double den  = 0.0;
    for( int i=1; i <= mx; i++ )
    {
        norm += fabs( qex.get(i,1) - qbc.get(i,1) );
        den  += fabs( qex.get(i,1) );
    }
    printf("\nStuff I've computed:\n");
    printf("  # of unkowns: %d;  Boundary condition order: %d\n", mx, mbc_order );
    printf("  The relative error with exact solution: %2.3e\n", norm / den );

//  for( int i=1; i <= mx; i++ )
//  {
//      printf("  q(%2d) = %2.3f\n", qbc.get(i,1) );
//  }

    return result;

}
