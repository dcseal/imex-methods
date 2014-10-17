#ifndef _REQ_PARAMS_H_
#define _REQ_PARAMS_H_
class req_params
{

    public:

        // starting and final time
        double tstart;
        double tfinal;

        // CFL numbers desired and acceptable
        double CFL_DES;
        double CFL_MAX;

        // number of equations/ghost cells
        int meqn;
        int maux;
        int mbc;
       
        // domain 
        double xlow;
        double xhigh;

        int mx;

        // spatial order
        int space_order;
       
        // type of ghost cells desired 
        int mbc_type_left;
        int mbc_type_right;

        double dtstart;        
        double dtmax;

        int nframes;

};

// init_params is a user defined function.  The underlying object here can be
// donwcast in order to recover user defined parameters
req_params* init_params();
void close_params( req_params* p );

#endif
