#ifndef _PARAMS_H_
#define _PARAMS_H_
#include "req_params.h"
class params : public req_params
{

    public:

        params();
        ~params();

        double u; // advection velocity

};
#endif
