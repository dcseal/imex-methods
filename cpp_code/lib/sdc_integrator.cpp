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
void IntegrateResidual( const dTensorBC3& Res, dTensorBC3& N);  

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
//
// This routine could be much more efficient if taylored to work with the
// butcher tableau rather than brute forcing its way through the tableau.
//
int sdc_integrator( req_state* top_state )
{

    dTensorBC2* q   = top_state->get_q();
    dTensorBC2* aux = top_state->get_aux();

    const double dt = top_state->get_dt();

    const int mx   = q->get_size(1);
    const int meqn = q->get_size(2);
    const int mbc  = q->get_mbc();

    const int num_intervals = 3;
    const int num_points    = num_intervals+1;

    // Itegrated Residual, right hand side function values and 
    // intermediate state values:
    dTensorBC3 N  ( mx, meqn, num_intervals, mbc );
    dTensorBC3 kE ( mx, meqn, num_points, mbc );
    dTensorBC3 kI ( mx, meqn, num_points, mbc );
    dTensorBC3 Q  ( mx, meqn, num_points, mbc );

    const double dt_small = dt / num_intervals;

    // Temporary storage
    dTensorBC2 qtmp( mx, meqn, mbc );
    for( int i=1-mbc; i <= mx+mbc; i++ )
    for( int m=1; m <= meqn; m++ )
    {
        Q.set(i, m, 1, q->get(i,m) );
        qtmp.set(i, m, q->get(i,m) );
    }

    // evaluate the right hand side:
    dTensorBC2 kItmp(mx, meqn, mbc );
    dTensorBC2 kEtmp(mx, meqn, mbc );
    ConstructL_Implicit  ( top_state->get_node(), &qtmp, aux, &kItmp );
    ConstructL           ( top_state->get_node(), &qtmp, aux, &kEtmp );
    for( int i=1-mbc; i <= mx+mbc; i++ )
    for( int m=1; m <= meqn; m++ )
    {
        kE.set(i,m,1, kEtmp.get(i,m) );
        kI.set(i,m,1, kItmp.get(i,m) );
    }

//  #ifdef IMPLICIT
//  // Used to store the rhs for the implicit solve
//  dTensorBC2 rhs( mx, meqn, mbc );
//  #endif

    // Take Forward/Backward  Euler steps with q:
    for( int s=1; s <= num_points-1; s++ )
    {

        // Forward Euler Steps
        for( int i=1-mbc; i <= mx+mbc; i++ )
        for( int m=1; m <= meqn; m++ )
        {
            const double frhs = kE.get(i, m, s) + kI.get(i,m,s);
            qtmp.set (i, m, qtmp.get(i,m) + dt_small*frhs );
            Q.set    (i, m, s+1, qtmp.get(i,m) );
        }

////    #ifdef IMPLICIT
////    // Backward Euler Step:
////    // yi = yi + dt A[i,i] * fI( yi );
////    if( s > 0 )
////    {


////        dTensor2* node = top_state->get_node();
////        rhs.CopyFrom( qtmp );

////        matrix_objs objs;

////        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
////           Compute the matrix and right-hand-side vector that define
////           the linear system, Ax = b.
////           - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

////        // TODO! //
////        const double eps = 1.0;
////        const double dx  = node->get(2,1)-node->get(1,1);

////        InitObjs( objs, mx, dx, g*eps*dt_small );

////        matrix_solve( objs, g*eps*dt_small, mbc, dx, rhs, qtmp );

////        /* 
////           Free work space.  All PETSc objects should be destroyed when they
////           are no longer needed.
////         */
////        CleanObjs    ( objs );


////    }
////    #endif

        // evaluate the right hand side:
        ConstructL_Implicit  ( top_state->get_node(), &qtmp, aux, &kItmp );
        ConstructL           ( top_state->get_node(), &qtmp, aux, &kEtmp );
        for( int i=1-mbc; i <= mx+mbc; i++ )
        for( int m=1; m <= meqn; m++ )
        {
            kE.set(i,m,s+1, kEtmp.get(i,m) );
            kI.set(i,m,s+1, kItmp.get(i,m) );
        }


    }

    const int num_corrections = 3;
    dTensorBC3 de( mx, meqn, num_points, mbc ); // TODO - this can be made smaller!
    de.setall(0.);
    for( int nc =1; nc <= num_corrections; nc++ )
    {

        // Step 1: Function Evaluation
        for( int n=1; n <= num_points; n++ ) // TODO - can start this at n=2
        {

            // Grab a q at stage level n
            for( int i=1-mbc; i <= mx+mbc; i++ )
            for( int m=1; m <= meqn; m++ )
            {
                double tmp = Q.get(i,m,n);
                qtmp.set(i,m,tmp);
            }

            // evaluate and save the right hand side:
            ConstructL_Implicit  ( top_state->get_node(), &qtmp, aux, &kItmp );
            ConstructL           ( top_state->get_node(), &qtmp, aux, &kEtmp );
            for( int i=1-mbc; i <= mx+mbc; i++ )
            for( int m=1; m <= meqn; m++ )
            {
                kE.set(i,m,n, kEtmp.get(i,m) );
                kI.set(i,m,n, kItmp.get(i,m) );
            }

        }

        // Step 2: Integrate the Residuatl
        IntegrateResidual( kE, N );


        // Step 3: March Delta forward through -
        // 
        // de^{k+1} = de^k - ( Q^{k+1} - Q^k ) + N^k + h ( f(Q^{k}+de^k) - f( Q^k ) )
        //

        for( int n=1; n <= num_intervals; n++ )
        {

            // First, evaluate f( Q^k + de^k ):
            for( int i=1-mbc; i <= mx+mbc; i++ )
            for( int m=1; m <= meqn; m++ )
            {
                double tmp = Q.get(i,m,n) + de.get(i,m,n);
                qtmp.set(i,m,tmp);
            }
            ConstructL( top_state->get_node(), &qtmp, aux, &kEtmp );


            // now we can march de forward
            for( int i=1-mbc; i <= mx+mbc; i++ )
            for( int m=1; m <= meqn; m++ )
            {
                double tmp = de.get(i,m,n) + 
                    0.5*dt*N.get(i,m,n) - ( Q.get(i,m,n+1)-Q.get(i,m,n) ) +
                    dt_small*( kEtmp.get(i,m) - kE.get(i,m,n) );
                de.set(i,m, n+1, tmp );
            }

        }

        // Step 4: Add the error into Q
        for( int i=1-mbc; i <= mx+mbc; i++ )
        for( int m=1; m <= meqn; m++ )
        for( int n=2; n <= num_points; n++ )
        {
            double tmp = de.get(i,m,n) + Q.get(i,m,n);
            Q.set(i,m,n,tmp);
        }

    }

    // Final value is what q should be set to:
    for( int i=1-mbc; i <= mx+mbc; i++ )
    for( int m=1; m <= meqn; m++ )
    {
        q->set(i,m, Q.get(i,m, num_points ) );
    }

    return 0;

}
