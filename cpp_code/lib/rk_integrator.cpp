#include <stdio.h>
#include <assert.h>
#include "defs.h"
#include "req_state.h"
#include "rk_coeffs.h"

#ifdef IMPLICIT
#include "petsc.h"
#include "matrix_solve.h"
#endif

void ConstructL( const dTensor2* node, dTensorBC2* q, dTensorBC2* aux, dTensorBC2* L );
void ConstructL_Implicit( const dTensor2* node, dTensorBC2* q, 
                                dTensorBC2* aux, dTensorBC2* L );

#ifdef IMPLICIT
int CleanObjs( matrix_objs & objs )
{
    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
     */
    PetscErrorCode ierr;
    ierr = VecDestroy(objs.x);CHKERRQ(ierr); 
    ierr = VecDestroy(objs.b);CHKERRQ(ierr); 

    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
     */
    ierr = MatDestroy(objs.A);CHKERRQ(ierr);
    ierr = KSPDestroy(objs.ksp);CHKERRQ(ierr);


}


int InitObjs( matrix_objs &objs, int mx, double dx, double as )
{

    Mat & A         = objs.A;
    KSP & ksp       = objs.ksp;
    PC  & pc        = objs.pc;

    Vec & x         = objs.x;
    Vec & b         = objs.b;

    PetscInt n = mx;
    PetscInt i, col[5];

    PetscErrorCode ierr;

    PetscScalar    value[5];
    /* 
       Create vectors.  Note that we form 1 vector from scratch and
       then duplicate as needed.
     */
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

    /* 
       Create matrix.  When using MatCreate(), the matrix format can
       be specified at runtime.

       Performance tuning note:  For problems of substantial size,
       preallocation of matrix memory is crucial for attaining good 
       performance. See the matrix chapter of the users manual for details.
     */
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);

    /* 
       Assemble matrix
     */
    const PetscScalar dxm2 = pow(dx,-2);
    value[0] = -1.0 * as * dxm2;
    value[1] =  2.0 * as * dxm2 + 1.0;
    value[2] = -1.0 * as * dxm2;
    for (i=1; i<n-1; i++) 
    {

        col[0] = i-1; col[1] = i; col[2] = i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

    }

    // apply the boundary conditions
    PetscScalar tmp = 0.;
    switch( 3 )
    {

        case 1:
        /* 1st order Dirichlet conditions on both sides */
        i = n - 1; col[0] = n - 2; col[1] = n - 1;
        //ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

        value[0] = 1.0;
        col[0]   = n-1;
        ierr = MatSetValues(A,1,&i,1,col,value,INSERT_VALUES);CHKERRQ(ierr);

        i = 0; col[0] = 0; col[1] = 1; value[0] = value[1]; value[1] = value[2];
        //ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

        value[0] = 1.0;
        col[0]   = 0;
        ierr = MatSetValues(A,1,&i,1,col,value,INSERT_VALUES);CHKERRQ(ierr);

        break;

        case 2:

        /* 2nd order dirichlet boundary conditions on both sides */
        i = n - 1; col[0] = n - 3; col[1] = n - 2; col[2] = n-1;

        value[0] = as * 1.0 * dxm2;
        value[1] = as * (-0.5) * dxm2;
        value[2] = as * (-7.0*0.5) * dxm2 + 1.0;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

        i = 0; col[0] = 0; col[1] = 1; col[2] = 2;
        value[2] = as * 1.0 * dxm2;
        value[1] = as * (-0.5) * dxm2;
        value[0] = as * (-7.0*0.5) * dxm2 + 1.0;

        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

        break;

        case 3: // periodic boundary conditions:

        i = n - 1; col[0] = n-2; col[1] = n - 1; col[2] = 0;
        tmp = value[1];
        value[0] = value[2]; 
        value[1] = value[0];
        value[2] = tmp;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

        //i = 0; col[0] = n-1; col[1] = 0; col[2] = 1;
        i = 0; col[1] = n-1; col[2] = 0; col[0] = 1;
        tmp = value[0];
        value[0] = value[1]; 
        value[1] = value[2];
        value[2] = tmp;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        break;

        default:
        printf("  bad\n");
        exit(1);

    }

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       Create the linear solver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Create linear solver context
     */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

    /* 
       Set operators. Here the matrix that defines the linear system
       also serves as the preconditioning matrix.
     */
    ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

    /* 
       Set linear solver defaults for this problem (optional).
       - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
       - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
     */
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

    /* 
       Set runtime options, e.g.,
       -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
       These options will override those specified above as long as
       KSPSetFromOptions() is called _after_ any other customization
       routines.
     */
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

 


}
#endif

// This function steps q from qn to qnp1.
// This function doesn't attempt to save q by any means.
// The choice of the RK coefficients is provided in rk_coeffs.h.
//
// This routine could be much more efficient if taylored to work with the
// butcher tableau rather than brute forcing its way through the tableau.
//
int rk_integrator( req_state* top_state, 
                        dTensorBC2** kE, dTensorBC2** kI, dTensorBC2& yi )
{

    dTensorBC2* q   = top_state->get_q();
    dTensorBC2* aux = top_state->get_aux();

    const double dt = top_state->get_dt();

    const int mx   = q->get_size(1);
    const int meqn = q->get_size(2);
    const int mbc  = q->get_mbc();

    dTensorBC2 rhs( mx, meqn, mbc );

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
                tmp += ( dt * AI[i][j] ) * kI[j]->get(ix,m);

            }
            yi.set( ix, m, tmp );

        }

        #ifdef IMPLICIT
        // implicit solve needs to happen here!!
        // yi = yi + dt A[i,i] * fI( yi );
        if( i > 0 )
        {


            dTensor2* node = top_state->get_node();
            rhs.CopyFrom( &yi );

            matrix_objs objs;

            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               Compute the matrix and right-hand-side vector that define
               the linear system, Ax = b.
               - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


            const double eps = 1.0;
            const double dx  = node->get(2,1)-node->get(1,1);

            InitObjs( objs, mx, dx, g*eps*dt );

            matrix_solve( objs, g*eps*dt, mbc, dx, rhs, yi );

            /* 
               Free work space.  All PETSc objects should be destroyed when they
               are no longer needed.
             */
            CleanObjs    ( objs );


        }
        #endif

        // evaluate the right hand side:
        ConstructL_Implicit  ( top_state->get_node(), &yi, aux, kI[i] );
        ConstructL           ( top_state->get_node(), &yi, aux, kE[i] );

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
            tmp+= dt * bI[n] * kI[n]->get( ix, m );
        }

        q->set( ix, m, tmp );

    }

    return 0;

/*
 * Depracated: Classical RK4 Integrator
 *
 */

/*
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
