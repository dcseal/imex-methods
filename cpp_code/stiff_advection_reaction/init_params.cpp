#include "params.h"
#include "math.h"

// User supplied function (really should be reading this from a text file!)
params::params()
{

    params* p = this;

    // global paramters
    p->tstart           = 0.0;
    p->tfinal           = 0.002;

    p->nframes          = 1;


    p->meqn             = 1;
    p->maux             = 1;
    p->mbc              = 2;
    p->mx               = 400;

    p->xlow             = 0.0;
    p->xhigh            = 1.0;

    const double dx = ( p->xhigh - p->xlow ) / p->mx;
    p->dtstart          = 0.5;
    p->dtmax            = 1.0e10;

    p->CFL_DES          = 0.9;
    p->CFL_MAX          = 1.0;

    p->space_order      = 4;

    p->mbc_type_left    = 0;
    p->mbc_type_right   = 0;

    // user specified parameters
    p->u                = 1.0;

}

params::~params()
{
    // the 'do nothing' destroyer!
    ;
}

// users need to do the proper casting here to delete extra user supplied
// parameters here
req_params* init_params()
{
    return (req_params*) new params();
}

void close_params( req_params* p )
{

    p = (params*) p;
    delete p;

}
void del_params(params* p)
{
    delete p;
}
