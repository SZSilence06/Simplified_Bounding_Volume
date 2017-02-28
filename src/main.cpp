
/*
 * Tool for generating simplified bounding volumes.
 * Author : SZ_Silence06
 * Date : Feb 24, 2017
 * Usage : sbvgen -s source_mesh_path [-d [output_directory]] [-e [max_distance]] [-options]
 * Possible options :
 *                     -t      Generate temp results.
 *                     -v      Display version information.
 */

#include <iostream>
#include <string.h>
#include <string>
#include "wkylib/CmdLine.h"

std::string g_inputMeshPath = "";
std::string g_outputPath = "";
double g_maxDistance = -1;
bool g_genTempResult = false;

using namespace WKYLIB;

void displayHelp()
{
    std::cout << "sbvgen tool for generating simplified bounding volumes." << std::endl
              << "Usage : sbvgen -s source_mesh_path [-d [output_directory]] [-e [max_distance]] [-options]" << std::endl
              << "Possible options :" << std::endl
              << "                    -t      Generate temp results." << std::endl
              << "                    -v      Display version information." <<std::endl;
}

void displayVersion()
{
    std::cout << "sbvgen version 1.00." << std::endl
              << "Author : SZ_Silence06" << std::endl
              << "Date : Feb 24, 2017"  << std::endl
              << "Source code available : https://github.com/SZSilence06/Simplified_Bounding_Volume" << std::endl
              << "Any problem please contact : SZ_Silence06@foxmail.com" << std::endl;
}

void parseCmdLines(int argc, char**argv)
{
    if(argc <= 1)
    {
        std::cout << "sbvgen version 1.00.\n"
                  << "Type -h for help."
                  << std::endl;
        exit(0);
    }

    CmdLine cmdParser(argc, argv);

    cmdParser.addParamDef("-s", CmdLine::CmdParamType::STRING);
    cmdParser.addParamDef("-d", CmdLine::CmdParamType::STRING);
    cmdParser.addParamDef("-e", CmdLine::CmdParamType::DOUBLE);
    cmdParser.addParamDef("-t", CmdLine::CmdParamType::BOOL);
    cmdParser.addParamDef("-h", CmdLine::CmdParamType::BOOL);
    cmdParser.addParamDef("-v", CmdLine::CmdParamType::BOOL);

    if(cmdParser.parse() == false)
    {
        std::cout << "Error occured while parsing input command!" << std::endl;
        switch(cmdParser.getErrorCode())
        {
        case CmdLine::ErrorCode::INVALID_PARAM_NAME:
            std::cout <<"Invalid parameter name: " << cmdParser.getErrorArg() << std::endl;
            break;
        case CmdLine::ErrorCode::INVALID_PARAM_VALUE:
            std::cout <<"Invalid parameter value: " << cmdParser.getErrorArg() << std::endl;
            break;
        case CmdLine::ErrorCode::NO_PARAM_VALUE:
            std::cout <<"No parameter value for parameter: " << cmdParser.getErrorArg() << std::endl;
            break;
        case CmdLine::ErrorCode::VALUE_OVERFLOW:
            std::cout <<"Paramter value overflows: " << cmdParser.getErrorArg() << std::endl;
            break;
        default:
            std::cout <<"Unknown error." << std::endl;
        }
        std::cout << "Type -h for help." << std::endl;
        exit(0);
    }

    bool needInputMesh = true;
    if(cmdParser.hasParam("-h"))
    {
        displayHelp();

        //now the user has requested to display help, so we dont't require user to input the mesh path.
        needInputMesh = false;
    }

    if(cmdParser.hasParam("-v"))
    {
        displayVersion();

        //now the user has requested to display version information, so we dont't require user to input the mesh path.
        needInputMesh = false;
    }

    if(cmdParser.hasParam("-s") == false && needInputMesh)
    {
        std::cout << "Cannot find '-s' parameter. Did you forget to type it before inputing the mesh file path?" << std::endl
                  << "Type -h for help." << std::endl;
        exit(0);
    }

    cmdParser.getString("-s", g_inputMeshPath);
    cmdParser.getString("-d", g_outputPath);
    cmdParser.getDouble("-e", g_maxDistance);
    cmdParser.getBool("-t", g_genTempResult);
}

int main(int argc, char**argv)
{
    parseCmdLines(argc, argv);

    return 0;
}
