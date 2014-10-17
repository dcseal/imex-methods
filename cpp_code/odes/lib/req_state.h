#ifndef __REQ_STATE_H_
#define __REQ_STATE_H_

#include "defs.h"
#include "req_params.h"
#include "user_funcs.h"

class req_state
{

    public:

        req_state( req_params* p );
        ~req_state( );

        // might as well declare these inline to speed up the code
        dTensorBC2* get_q()   const        { return q;       }
        dTensorBC2* get_aux() const        { return aux;     }
        dTensor2* get_node()  const        { return node;    }

        double get_dtstart()  const        { return dtstart; }
        double get_t ( ) const             { return this->t;  }
        double get_dt( ) const             { return this->dt; }

        void set_dtstart( double dt )      { dtstart = dt;   }
        void set_t( double t )             { this->t = t;    }
        void set_dt( double dt )           { this->dt = dt;  }

        // NOTE: THIS WILL BE TYPECAST DOWN IF USER HAS ANY LOCAL PARAMS THEY
        // WISH TO USE
        req_params* get_params() const { return p;      }

    protected:
        virtual void set_grid();

    private:

        
        // stuff declared with a new operator
        dTensorBC2* q;
        dTensorBC2* aux;
        req_params* p;
        dTensor2* node;

        // stuff we don't have to destroy
        double t;
        double dtstart;
        double dt;

};

// user defined routines here for allocating and deleting state variables
req_state* init_state(req_params* p);
void close_state( req_state* top_state );

#endif
