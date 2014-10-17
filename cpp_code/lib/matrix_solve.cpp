#include <stdio.h>
#include <stdlib.h>
#include "tensors.h"
#include "petscksp.h"
#include "matrix_solve.h"

void ExtractVec( Vec v, dTensorBC2& q )
{
    /* Extract the PETSC vec object into my object */

    const int mx   = q.get_size(1);
    const int meqn = q.get_size(2);

    double mq = 0.;
    for( int i=1; i <= mx;   i++ )
    for( int m=1; m <= meqn; m++ )
    {

        PetscInt    ix[1], ni=1;
        PetscScalar qi[1];

        //TODO - CHANGE THIS INDEXING!
        ix[0] = i-1;
        VecGetValues(v,ni,ix,qi);

        double tmp = qi[0];
        q.set(i, m, tmp );
    }

}

int VecInit( matrix_objs & objs, PetscInt n )
{

    Vec & x         = objs.x;
    Vec & b         = objs.b;

    /* 
       Create vectors.  Note that we form 1 vector from scratch and
       then duplicate as needed.
     */
    PetscErrorCode ierr;
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

}

int MatrixInit(matrix_objs& objs, int mx, double dx, double as )
{

    //////////////////////////////////////////////////////////////////////////
    //////////// local variables / references ////////////
    Mat & A         = objs.A;
    KSP & ksp       = objs.ksp;
    PC  & pc        = objs.pc;


    PetscErrorCode ierr;
    PetscInt n = mx;
    PetscInt       i, col[3];
    PetscScalar    value[5];
    //////////////////////////////////////////////////////////////////////////

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
        printf("  bad choice for boundary conditions on matrix.\n");
        exit(1);

    }

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

}


int matrix_solve( matrix_objs& objs, double as, const int mbc_order,
    double dx, const dTensorBC2& f, dTensorBC2& q )
{


    Mat & A         = objs.A;
    KSP & ksp       = objs.ksp;
    PC  & pc        = objs.pc;

    Vec & x         = objs.x;
    Vec & b         = objs.b;

    PetscReal      norm;         /* norm of solution error */
    PetscErrorCode ierr;

    const int mx   = q.get_size(1);
    const int meqn = q.get_size(2);
    const int mbc  = q.get_mbc();

    PetscInt       n = mx;

    PetscInt       i, col[3],its;

    PetscMPIInt    size;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(1,"This is a uniprocessor example only!");
    ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);

    PetscInt ix[1];
    for (i=1; i<=mx; i++) 
    {
        // TODO - fix the indexing here ... //
        ierr = VecSetValue(b, i-1, f.get(i,1), INSERT_VALUES); CHKERRQ(ierr);
    }

  
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
    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       Check solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const bool verbose = false;
    if( verbose )
    {
        ierr = KSPGetIterationNumber(ksp,&its);     CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);CHKERRQ(ierr);
    }

    ExtractVec( x, q );

    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
     */
//  CleanObjs    ( objs );

    return 0;

}
