#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main( int argc, char** args )
{

    char command_str[1024];
    //const char* get_outputdir();
    char outputdir[] = "/home/seal4/im/cpp_code/scripts/output1/output2/output3";
    
    snprintf( command_str, 1024,
        "./create_output_dir.py -o %s\n", outputdir );

    return system( command_str );

}

