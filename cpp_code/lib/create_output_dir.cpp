#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <stdlib.h>
#include <stdio.h>

#include "req_params.h"

using namespace std;

int create_output_dir( string outputdir, req_params* p )

{

    int nframes, meqn, melems;
    double xlow, xhigh, dx;

    nframes = p->nframes;
    meqn    = p->meqn;
    melems  = p->mx;
    xlow    = p->xlow;
    xhigh   = p->xhigh;
    dx      = (p->xhigh - p->xlow) / melems;

    char command_str[1024];
    snprintf( command_str, 1024,
        "$FIN_VOL_DIR/scripts/create_output_dir.py -o %s\n", outputdir.c_str() );

    // if this doesn't return 0, something bad happened //
    int sys_result = system( command_str );

    // Output meqn and nout for plotting purposes
    string qhelp;

    qhelp=outputdir+"/qhelp.dat";
    ofstream out_file(qhelp.c_str(), ios::out);

    out_file << setprecision(16);
    out_file << nframes << endl << meqn << endl << melems << endl;
    out_file << setw(24) << scientific << xlow << endl << setw(24)
        << scientific << xhigh << endl << setw(24) 
        << scientific << dx << endl;
    out_file.close();

}
