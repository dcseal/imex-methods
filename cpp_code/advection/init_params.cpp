#include "params.h"

// User supplied function (really should be reading this from a text file!)
params::params()
{

    params* p = this;

    // global paramters
    tstart           = 0.0;
    tfinal           = 1.0;

    nframes          = 1;


    CFL_DES          = 0.99;
    CFL_MAX          = 1.0;

    meqn             = 1;
    maux             = 0;
    mbc              = 2;
    mx               = 400;

    xlow             = 0.0;
    xhigh            = 1.0;

    dtstart          = 2.0 * (p->xhigh - p->xlow) / p->mx;
    dtmax            = 1.0e10;

    space_order      = 4;

    mbc_type_left    = 0;
    mbc_type_right   = 0;

    // user specified parameters
    u                = 1.0;

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
