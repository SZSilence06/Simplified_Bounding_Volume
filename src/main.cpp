
/*
 * Tool for generating simplified bounding volumes.
 * Author : SZ_Silence06
 * Date : Feb 24, 2017
 * Usage : sbvgen [source mesh] [-d [output directory]] [-e [max distance]] [-options]
 * Possible options:
 *                  -t      Generate temp results.
 */

#include <iostream>
#include <string.h>
#include <string>

std::string g_inputMeshPath;
std::string g_outputPath;

void parseCmdLines(int argc, char**argv)
{
    char* cmd = nullptr;
    char* param = nullptr;
    if(argc <= 1)
    {
        std::cout << "sbvGen version 1.00.\n"
                  << "Type -h for help."
                  << std::endl;
        exit(0);
    }
    g_inputMeshPath = std::string(argv[1]);
}

int main(int argc, char**argv)
{
    parseCmdLines(argc, argv);



    return 0;
}
