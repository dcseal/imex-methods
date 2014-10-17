#include "math.h"
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "tensors.h"
#include "defs.h"
using namespace std;

void output( const dTensorBC2& q, double t, int nframe)
{

    const int melems  = q.get_size(1);
    const int meqn    = q.get_size(2);
    const int mbc     = q.get_mbc();

    string outputdir("./output/");

    // Open file -- q
    ostringstream fname1;
    fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream q_file(fname1.str().c_str(), ios::out );

    q_file << setprecision(16);
    q_file << setw(24) << scientific << t << endl;

    // Output each coefficient
        for (int m=1; m<=meqn; m++)
            for (int i=1; i<=melems; i++)      
            {
                double tmp = q.get(i,m);

                if (fabs(tmp) < 1.0e-99) {tmp = 0.0;}
                q_file << setw(24) << scientific << tmp << endl;
            }
    q_file.close();

    // Open file -- aux
//  ostringstream fname2;
//  fname2 << outputdir << "/" << "a" << setfill('0') 
//      << setw(4) << nframe << ".dat";
//  ofstream aux_file(fname2.str().c_str(), ios::out );

//  aux_file << setprecision(16);
//  aux_file << setw(24) << scientific << t << endl;

//  // Output each coefficient
//  for (int k=1; k<=method1; k++)
//      for (int m=1; m<=maux; m++)
//          for (int i=1; i<=melems; i++)      
//          {
//              double tmp = aux.get(i,m,k);

//              if (fabs(tmp) < 1.0e-99) {tmp = 0.0;}
//              aux_file << setw(24) << scientific << tmp << endl;
//          }
//  aux_file.close();
}
