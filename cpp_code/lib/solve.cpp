#include "defs.h"
#include "CONSTANTS.h"
#include "req_state.h"

#include "rk_coeffs.h"
//#define MAX_STAGES ( 6 )
int rk_integrator( req_state* top_state, 
                        dTensorBC2** kE, dTensorBC2** kI, dTensorBC2& yi );

int sdc_integrator( req_state* top_state );


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

    const double CFL_MAX = top_state->get_CFL_MAX();
    const double CFL_DES = top_state->get_CFL_DES();

    dTensorBC2 qold(mx,meqn,mbc);

    ///////////////////////////////////////////////////////////////////////////
    // temporary storage ( this should be called from above! )
    dTensorBC2 q_tmp( mx, meqn, mbc );

    dTensorBC2** kE = new dTensorBC2*[ (s) * sizeof( dTensorBC2* ) ];
    dTensorBC2** kI = new dTensorBC2*[ (s) * sizeof( dTensorBC2* ) ];
    for( int j=0; j < s; j++ )
    {   
        kE[j] = new dTensorBC2( mx, meqn, mbc );
        kI[j] = new dTensorBC2( mx, meqn, mbc ); 
    }
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
     
            top_state->set_dt( dt );
            //rk_integrator( top_state, kE, kI, q_tmp );
            sdc_integrator( top_state );

// TODO - DO WE CARE ABOUT DOING ANYTHING FANCY WITH CFL NUMBERS HERE?  WHAT
// ABOUT REJECTING STEPS?
//          cfl = get_cfl( smax, q );
//const double umax = p->u;
const double umax = 1.0;
double cfl = umax * dt / dx;

            // chose a new time step:
            if( cfl > 0.0)
            {
                dt = min( dtv[2], dt*CFL_DES/cfl );
                printf("  dt = %2.3e; cfl = %2.3f; t = %2.3e\n", dt, cfl, t );
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
                m_accept = 0;
            }

        }

    }


    // save the time step for next call to this function
    top_state->set_dtstart( dt );


    ///////////////////////////////////////////////////////////////////////////
    // cleanup all the temporary storage used
    for( int j=0; j < s; j++ )
    { 
        delete kE[j]; 
        delete kI[j];
    }
    delete[] kE;
    delete[] kI;
    ///////////////////////////////////////////////////////////////////////////


}
