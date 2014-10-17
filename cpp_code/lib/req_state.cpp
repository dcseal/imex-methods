#include "req_state.h"

// note: most of these function are simply declared inline
req_state::req_state( req_params* p )
{

    this->q    = new dTensorBC2( p->mx, p->meqn, p->mbc );
    this->aux  = new dTensorBC2( p->mx, p->maux, p->mbc );
    this->node = new dTensor2  ( p->mx + 1, 1 );
    this->p    = p;

    set_grid();    

}


req_state::~req_state( )
{

    delete this->q;
    delete this->aux;
    delete this->node;

}

void req_state::set_grid()
{


    // create the grid
    const int mx   = p->mx;
    const int meqn = p->meqn;
    const int mbc  = p->mbc;

    const double xlow = p->xlow;
    const double dx   = (p->xhigh - p->xlow) / mx;
    for( int n=1; n <= mx+1; n++ )
    {
        this->node->set(n, 1, xlow + (n-1) * dx );
    }

}
