#include "defs.h"
#include "CONSTANTS.h"
#include "req_state.h"

void rk_integrator( req_state* top_state, 
                        dTensorBC2** kE, dTensorBC2& yi );

// This routine is responsible for solving the problem from two input times
// tstart to tend.
void solve( req_state* top_state, double Delta_t )
{



    dTensorBC2* q   = top_state->get_q();

    double t        = top_state->get_t();
    double dt       = top_state->get_dtstart();
    double dtv[]    = { 0.0, dt, 1.0 };

    const int mx    = q->get_size(1);
    const int meqn  = q->get_size(2);
    const int mbc   = q->get_mbc();

    const dTensor2* node = top_state->get_node();

    const double dx = node->get(2,1) - node->get(1,1);

//  const double CFL_MAX = this->p->CFL_MAX;
//  const double CFL_DES = this->p->CFL_DES;
    const double CFL_MAX = 0.95;
    const double CFL_DES = 0.9;

    dTensorBC2 qold(mx,meqn,mbc);

    ///////////////////////////////////////////////////////////////////////////
    // temporary storage ( this should be called from above! )
    dTensorBC2 q_tmp( mx, meqn, mbc );
const int s = 6;
    dTensorBC2** kE = new dTensorBC2*[ (s) * sizeof( dTensorBC2* ) ];
    for( int j=0; j < s; j++ )
    {   kE[j] = new dTensorBC2( mx, meqn, mbc ); }
    ///////////////////////////////////////////////////////////////////////////



    // main time integration loop
    const double tend = t + Delta_t;
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
     
            //rk_integrator( p->xlow, dx, dt, kE, q_tmp, q );
            top_state->set_dt( dt );
            rk_integrator( top_state, kE, q_tmp );

// TODO - DO WE CARE ABOUT DOING ANYTHING FANCY WITH CFL NUMBERS HERE?  WHAT
// ABOUT REJECTING STEPS?
//          cfl = get_cfl( smax, q );
//const double umax = p->u;
const double umax = 1.0;

//          double cfl = p->u * dt / dx;
double cfl = umax * dt / dx;

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
    top_state->set_dtstart( dt );


    ///////////////////////////////////////////////////////////////////////////
    // cleanup all the temporary storage used
    for( int j=0; j < s; j++ )
    { delete kE[j]; }
    delete[] kE;
    ///////////////////////////////////////////////////////////////////////////


}
