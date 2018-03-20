
/*
 * Tool for generating simplified bounding volumes.
 * Author : SZ_Silence06
 * Date : Feb 24, 2017
 * Usage : sbvgen -s source_mesh_path [-d [output_directory]] [-e [max_distance]] [-r [sample radius]] [-a alpha_param_value] [-options]
 * Possible options :
 *                     -t      Generate intermediate results.
 *                     -v      Display version information.
 *                     -f      Visualize field instead of simplifying the geometry.
 */

#include <iostream>
#include <string.h>
#include <string>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>
#include <wkylib/CmdLine.h>
#include <wkylib/mesh/IO.h>
#include <wkylib/debug_util.h>
#include "alg/Simplifier.h"
#include "alg/ShellGenerator.h"
#include "alg/Shell.h"
#include "alg/Logger.h"

std::string g_inputMeshPath = "";
std::string g_outputPath = "";
double g_maxDistance = -1;
double g_sampleRadius = -1;
bool g_genTempResult = false;
bool g_visualizeField = false;
double g_alpha = std::numeric_limits<double>::max();

using namespace WKYLIB;

void displayHelp()
{
    std::cout << "sbvgen tool for generating simplified bounding volumes." << std::endl
              << "Usage : sbvgen -s source_mesh_path [-d [output_directory]] [-e [max_distance]] "
                 "[-r [sample radius]] [-options]" << std::endl
              << "Possible options :" << std::endl
              << "                    -t      Generate temp results." << std::endl
              << "                    -v      Display version information." <<std::endl
              << "                    -f      Visualize field instead of simplifying the geometry." << std::endl;
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
    cmdParser.addParamDef("-r", CmdLine::CmdParamType::DOUBLE);
    cmdParser.addParamDef("-a", CmdLine::CmdParamType::DOUBLE);
    cmdParser.addParamDef("-t", CmdLine::CmdParamType::BOOL);
    cmdParser.addParamDef("-h", CmdLine::CmdParamType::BOOL);
    cmdParser.addParamDef("-v", CmdLine::CmdParamType::BOOL);
    cmdParser.addParamDef("-f", CmdLine::CmdParamType::BOOL);

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

    if(needInputMesh == false)
    {
        exit(0);
    }

    if(cmdParser.hasParam("-s") == false)
    {
        std::cout << "Cannot find '-s' parameter. Did you forget to type it before inputing the mesh file path?" << std::endl
                  << "Type -h for help." << std::endl;
        exit(0);
    }

    cmdParser.getString("-s", g_inputMeshPath);
    cmdParser.getString("-d", g_outputPath);
    cmdParser.getDouble("-e", g_maxDistance);
    cmdParser.getDouble("-r", g_sampleRadius);
    cmdParser.getDouble("-a", g_alpha);
    cmdParser.getBool("-t", g_genTempResult);
    cmdParser.getBool("-f", g_visualizeField);
}

void genDefaultParams()
{
    if(g_alpha == std::numeric_limits<double>::max())
    {
        g_alpha = 0.2;
    }

    if(g_outputPath == "")
    {
        //no output directory, so we auto generate one.
        size_t pos = g_inputMeshPath.find_last_of('.');
        if(pos == std::string::npos)
        {
            g_outputPath = g_inputMeshPath;
        }
        else
        {
            g_outputPath = g_inputMeshPath.substr(0, pos);
        }
        std::stringstream ss;
        ss << g_outputPath << "_-e_ " << g_maxDistance << "_-r_" << g_sampleRadius << "_-a_" << g_alpha;
        g_outputPath = ss.str();
    }

    if(boost::filesystem::exists(g_outputPath) == false
            || boost::filesystem::is_directory(g_outputPath) == false)
    {
        //output directory don't exist, so create it
        boost::filesystem::create_directory(g_outputPath);
    }
}

void generateShells(const SBV::Mesh& mesh, SBV::Shell& shell)
{
    WKYLIB::DebugTimer timerGenerateShell("Generate Shell");
    timerGenerateShell.start();

    SBV::ShellGenerator generator(mesh.vertices, mesh.triangles, g_outputPath);
    generator.generate(g_maxDistance, g_sampleRadius, shell, g_visualizeField);

    //WKYLIB::Mesh::readPoints(g_outputPath + "/inner_shell.vtk", shell.mInnerShell);
    //WKYLIB::Mesh::readPoints(g_outputPath + "/outer_shell.vtk", shell.mOuterShell);
    //shell.buildKdTree();

    if(g_genTempResult)
    {
        WKYLIB::Mesh::writePoints(g_outputPath + "/inner_shell.vtk", shell.mInnerShell);
        WKYLIB::Mesh::writePoints(g_outputPath + "/outer_shell.vtk", shell.mOuterShell);
    }

    timerGenerateShell.end();
    SBV::Logger& logger = SBV::Logger::getInstance();
    logger.setFile(g_outputPath + "/log.txt");
    logger.log("Generate Shell : " + std::to_string(timerGenerateShell.getTime()) + " ms.");
}

int main(int argc, char**argv)
{
    parseCmdLines(argc, argv);

    SBV::Mesh mesh;
    if(jtf::mesh::load_obj(g_inputMeshPath.c_str(), mesh.triangles, mesh.vertices))
    {
        std::cout << "Fail to load mesh: " + g_inputMeshPath << std::endl;
        return 0;
    }

    genDefaultParams();

    SBV::Shell shell;
    generateShells(mesh, shell);

    SBV::Simplifier simplifier(shell);
    simplifier.setOutputDirectory(g_outputPath);
    simplifier.setAlpha(g_alpha);
    simplifier.setGenTempResult(g_genTempResult);

    simplifier.simplify();

    return 0;
}
