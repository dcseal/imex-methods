#include "req_state.h"
#include "params.h"

// users wishing to add extra parameters should replace these functions
// with their own desired functions.
req_state* init_state(req_params* p)
{
    return new req_state( p );
}

void close_state( req_state* st )
{
    delete st;
}
