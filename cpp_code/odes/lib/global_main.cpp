#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "req_state.h"
#include "req_params.h"
#include "user_funcs.h"

void output( const dTensorBC2& q, double t, int framenum );
void solve( req_state* top_state, double Delta_t );


using namespace std;

int global_main( )
{

    // Print a pretty welcome screen //
    printf("///////////////////////////////////////");
    printf("///////////////////////////////////////\n");

    printf("//                                     ");
    printf("                                     //\n");

    printf("// Starting the Program:               ");
    printf("                                     //\n");

    printf("//                                     ");
    printf("                                     //\n");

    printf("// Maybe one day I'll put the parameter");
    printf("s here                               //\n");

    printf("///////////////////////////////////////");
    printf("///////////////////////////////////////\n");


    // set all of the user defined parameters:
    req_params* p        = init_params();
    req_state* top_state = init_state( p );


    // initial condition for q and aux
    const int mx  = p->mx;
    const int mbc = p->mbc;
    void integrate_on_cells( int sorder, 
        const dTensor2& node, dTensorBC2& q,  
        void (*Func)( const dTensor2& xpts, dTensor2& qpts ) );
    integrate_on_cells( 4, *top_state->get_node(), 
                            *top_state->get_q(), Qinit );

    if( p->maux > 0 )
    {
        integrate_on_cells( 4, *top_state->get_node(), 
                            *top_state->get_aux(), AuxInit );
    }

    // helper file for the output directory
    string outputdir;
    outputdir = "./output/";

    void create_output_dir( string outputdir,
        int nframes, int meqn, int melems, double xlow, double xhigh, double dx);
    create_output_dir( outputdir, p->nframes, p->meqn, p->mx, p->xlow, 
        p->xhigh, (p->xhigh - p->xlow) / p->mx );



    // perform a solve for each frame number
    const int nframes = p->nframes;
    output( *top_state->get_q(), 0.0, 0);      // output the initial conditions

    double t  = p->tstart;
    double Delta_t = (p->tfinal - p->tstart) / nframes;
    for( int n=1; n <= p->nframes; n++ )
    {
        solve( top_state, Delta_t );
        output( *top_state->get_q(), top_state->get_t(), n ); // frome number == n

        printf("frame number: %3d ; t = %2.3e ; \n", n, t );
    }

    close_state( top_state );
    close_params( p );

}
