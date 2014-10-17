#include <stdio.h>
#include <stdlib.h>
#include "tensors.h"
#include "petscksp.h"

void ExtractVec( Vec v, dTensorBC2& q )
{
    /* Extract the PETSC vec object into my object */

    const int mx   = q.get_size(1);
    //const int meqn = q.get_size(2);

    for( int i=1; i <= mx; i++ )
    {

        PetscInt    ix[1], ni=1;
        PetscScalar qi[1];
        ix[0] = i-1;
        VecGetValues(v,ni,ix,qi);

        double tmp = qi[0];
        q.set(i, 1, tmp );
    }

}

int matrix_solve( const int mbc_order, const dTensorBC2& qex, 
    double dx, const dTensorBC2& f, dTensorBC2& q, int argc, char** args)
{


    Vec            x, b, u;      /* approx solution, RHS, exact solution */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    PetscReal      norm;         /* norm of solution error */
    PetscErrorCode ierr;

    const int mx   = q.get_size(1);
    const int meqn = q.get_size(2);
    const int mbc  = q.get_mbc();

    PetscInt       n = mx;

    PetscInt       i, col[3],its;

    PetscMPIInt    size;
    PetscScalar    neg_one = -1.0,one = 1.0,value[5];

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(1,"This is a uniprocessor example only!");
    ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       Compute the matrix and right-hand-side vector that define
       the linear system, Ax = b.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /* 
       Create vectors.  Note that we form 1 vector from scratch and
       then duplicate as needed.
     */
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

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
    value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
    for (i=1; i<n-1; i++) 
    {

        col[0] = i-1; col[1] = i; col[2] = i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

    }

    switch( mbc_order )
    {

        case 1:
        /* 1st order Dirichlet conditions on both sides */
        i = n - 1; col[0] = n - 2; col[1] = n - 1;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

        i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

        break;

        case 2:

        /* 2nd order dirichlet boundary conditions on both sides */
        i = n - 1; col[0] = n - 3; col[1] = n - 2; col[2] = n-1;
        value[0] = 1.0; value[1] = -0.5; value[2] = -7.0*0.5;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        i = 0; col[0] = 0; col[1] = 1; col[2] = 2;
        value[2] = 1.0; value[1] = -0.5; value[0] = -7.0*0.5;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

        break;

        default:
        printf("  bad\n");
        exit(1);

    }

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* 
       Set exact solution; then compute right-hand-side vector.
     */
//  ierr = VecSet(u,one);CHKERRQ(ierr);
//  ierr = MatMult(A,u,b);CHKERRQ(ierr);

    PetscInt ix[1];
    PetscScalar y[1];
    for (i=1; i<=mx; i++) 
    {
        //ierr = VecSetValues(b, 1, ix, y, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(b, i-1, f.get(i,1)*pow(dx,2), INSERT_VALUES); CHKERRQ(ierr);

        //ierr = VecSetValues(u, 1, ix, y, INSERT_VALUES); CHKERRQ(ierro);
        ierr = VecSetValue(u, i-1, qex.get(i,1), INSERT_VALUES); CHKERRQ(ierr);
    }



    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       Create the linear solver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Create linear solver context
     */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

    // This part is broken when using non-symmetric matrix!
    //KSPSetType(ksp, KSPCG);

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
    //ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    //ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    //ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,100000);CHKERRQ(ierr);

    /* 
       Set runtime options, e.g.,
       -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
       These options will override those specified above as long as
       KSPSetFromOptions() is called _after_ any other customization
       routines.
     */
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       Solve the linear system
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Solve linear system
     */
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr); 

    /* 
       View solver info; we could instead use the option -ksp_view to
       print this info to the screen at the conclusion of KSPSolve().
     */
    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    ExtractVec( x, q );

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       Check solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Check the error
     */
    ierr = VecAXPY(x,neg_one,u);                CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&norm);             CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);     CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %A, Iterations %D\n",
            norm,its);CHKERRQ(ierr);


    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
     */
    ierr = VecDestroy(x);CHKERRQ(ierr); ierr = VecDestroy(u);CHKERRQ(ierr);
    ierr = VecDestroy(b);CHKERRQ(ierr); ierr = MatDestroy(A);CHKERRQ(ierr);
    ierr = KSPDestroy(ksp);CHKERRQ(ierr);

    /*
       Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
       options are chosen (e.g., -log_summary).
     */

    return 0;

}
