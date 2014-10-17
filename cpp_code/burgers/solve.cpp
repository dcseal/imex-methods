#include "defs.h"
#include "app_solver.h"
#include "CONSTANTS.h"

void rk_integrator( const double xlow, const double dx, const double dt, 
                        dTensorBC2** kE, dTensorBC2& yi, dTensorBC2* q );

// This routine is responsible for solving the problem from two input times
// tstart to tend.
void app_solver::solve( double tstart, double tend )
{


    double t        = tstart;
    double dt       = this->dtstart;
    double dtv[]    = { 0.0, this->dtstart, this->p->dtmax };
    dTensorBC2* q   = this->q;

    const int mx    = q->get_size(1);
    const int meqn  = q->get_size(2);
    const int mbc   = q->get_mbc();

    const double dx = ( p->xhigh - p->xlow ) / mx;

    const double CFL_MAX = this->p->CFL_MAX;
    const double CFL_DES = this->p->CFL_DES;

    dTensorBC2 qold(mx,meqn,mbc);

    ///////////////////////////////////////////////////////////////////////////
    // temporary storage ( this should be called from above! )
// TODO -- THIS SHOULDN'T BE HARD CODED HERE! //
const int s = 6;
    dTensorBC2 q_tmp( mx, meqn, mbc );
    dTensorBC2** kE = new dTensorBC2*[ (s) * sizeof( dTensorBC2* ) ];
    for( int j=0; j < s; j++ )
    {   kE[j] = new dTensorBC2( mx, meqn, mbc ); }
    ///////////////////////////////////////////////////////////////////////////



    // main time integration loop
    while( t < tend )
    {

       
        // copy in case we need to reject this step
        qold.CopyFrom( q );

        int m_accept = 0;
        while( m_accept == 0 )
        {

            double told = t;
            if( told + dt > tend )
            { dt = tend-told; }
            t = t + dt; 
     
            rk_integrator( p->xlow, dx, dt, kE, q_tmp, q );

// TODO - DO WE CARE ABOUT DOING ANYTHING FANCY WITH CFL NUMBERS HERE?  WHAT
// ABOUT REJECTING STEPS?
//          cfl = get_cfl( smax, q );
const double umax = p->u;

            double cfl = p->u * dt / dx;

            // output time step information
            /// TODO ! //


            // chose a new time step:
            if( cfl > 0.0)
            {
                dt = min( dtv[2], dt*CFL_DES/cfl );
//              double dtmin = min(dt,dtmin);
//              double dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            if( cfl <= CFL_MAX ) // accept!
            { m_accept = 1; }
            else // reject!
            {
                t  = told;
                q->CopyFrom( &qold );
                printf("  rejecting step! dt = %2.5e\n", dt);
            }

        }

    }


    // save the time step for next call to this function
    this->dtstart = dt;


    ///////////////////////////////////////////////////////////////////////////
    // cleanup all the temporary storage used
    for( int j=0; j < s; j++ )
    { delete kE[j]; }
    delete[] kE;
    ///////////////////////////////////////////////////////////////////////////


}
