
/*
 * Tool for generating simplified bounding volumes.
 * Author : SZ_Silence06
 * Date : Feb 24, 2017
 * Usage :
 *     There are two categories of usages. One will generate shells, another build the simplified mesh with the generated shells.
 *
 *     (1) For generating shells, usage is:
 *
 *     sbvgen --shell -s source_mesh_path [-d output_directory] [-e max_distance] [-r sample_radius] [-options]
 *
 *     example : sbvgen --shell -s bunny_2.obj -e 0.05 -r 0.05
 *
 *     Posiible options:
 *
 *                     -f      Visualize field instead of simplifying the geometry.
 *
 *     (2) For simplification, usage is:
 *
 *     sbvgen --inner inner_shell_file_path --outer outer_shell_file_path [-a alpha_param_value] [-options]
 *
 *     example : sbvgen --inner inner_shell.vtk --outer outer_shell.vtk -t
 *
 *     Possible options :
 *                     -t      Generate intermediate results.
 *                     -v      Display version information.
 *
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
std::string g_innerSamplePath;
std::string g_outerSamplePath;
double g_maxDistance = -1;
double g_sampleRadius = -1;
bool g_genTempResult = false;
bool g_visualizeField = false;
bool g_genShell = false;
double g_alpha = std::numeric_limits<double>::max();

using namespace WKYLIB;

void displayHelp()
{
    std::cout << "Tool for generating simplified bounding volumes." << std::endl
              << "Author : SZ_Silence06" << std::endl
              << "Date : Feb 24, 2017" << std::endl
              << "Usage :" << std::endl
              << "     There are two categories of usages. One will generate shells, another build the simplified mesh with the generated shells." << std::endl
              << std::endl
              << "     (1) For generating shells, usage is:" << std::endl
              << std::endl
              << "     sbvgen --shell -s source_mesh_path [-d output_directory] [-e max_distance] [-r sample_radius] [-options]" << std::endl
              << std::endl
              << "     example : sbvgen --shell -s bunny_2.obj -e 0.05 -r 0.05" << std::endl
              << std::endl
              << "     Posiible options:" << std::endl
              << std::endl
              << "                     -f      Visualize field instead of simplifying the geometry." << std::endl
              << std::endl
              << "     (2) For simplification, usage is:" << std::endl
              << std::endl
              << "     sbvgen --inner inner_shell_file_path --outer outer_shell_file_path [-a alpha_param_value] [-options]" << std::endl
              << std::endl
              << "     example : sbvgen --inner inner_shell.vtk --outer outer_shell.vtk -t" << std::endl
              << std::endl
              << "     Possible options :" << std::endl
              << "                     -t      Generate intermediate results." << std::endl
              << "                     -v      Display version information." << std::endl;
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
    cmdParser.addParamDef("--shell", CmdLine::CmdParamType::BOOL);
    cmdParser.addParamDef("--inner", CmdLine::CmdParamType::STRING);
    cmdParser.addParamDef("--outer", CmdLine::CmdParamType::STRING);

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

    cmdParser.getBool("--shell", g_genShell);
    if(g_genShell) {
        if(cmdParser.hasParam("-s") == false)
        {
            std::cout << "Cannot find '-s' parameter. Did you forget to type it before inputing the mesh file path?" << std::endl
                      << "Type -h for help." << std::endl;
            exit(0);
        }
        cmdParser.getString("-s", g_inputMeshPath);
        cmdParser.getDouble("-e", g_maxDistance);
        cmdParser.getDouble("-r", g_sampleRadius);
        cmdParser.getBool("-f", g_visualizeField);
    }
    else {
        if(cmdParser.hasParam("--inner") == false)
        {
            std::cout << "Cannot find '--inner' parameter. Did you forget to type it before inputing the mesh file path?" << std::endl
                      << "Type -h for help." << std::endl;
            exit(0);
        }
        cmdParser.getString("--inner", g_innerSamplePath);

        if(cmdParser.hasParam("--inner") == false)
        {
            std::cout << "Cannot find '--outer' parameter. Did you forget to type it before inputing the mesh file path?" << std::endl
                      << "Type -h for help." << std::endl;
            exit(0);
        }
        cmdParser.getString("--outer", g_outerSamplePath);

        cmdParser.getDouble("-a", g_alpha);
        cmdParser.getBool("-t", g_genTempResult);
    }

    cmdParser.getString("-d", g_outputPath);      
}

void genDefaultParams()
{
    if(g_alpha == std::numeric_limits<double>::max())
    {
        g_alpha = 0.2;
    }

    if(g_outputPath == "")
    {
        if(g_genShell){
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
            ss << g_outputPath << "_-e_" << g_maxDistance << "_-r_" << g_sampleRadius;
            g_outputPath = ss.str();
        }
        else {
            g_outputPath = ".";
        }
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

    WKYLIB::Mesh::writePoints(g_outputPath + "/inner_shell.vtk", shell.mInnerShell);
    WKYLIB::Mesh::writePoints(g_outputPath + "/outer_shell.vtk", shell.mOuterShell);

    timerGenerateShell.end();
    SBV::Logger& logger = SBV::Logger::getInstance();
    logger.log("Generate Shell : " + std::to_string(timerGenerateShell.getTime()) + " ms.");
}

void readShell(SBV::Shell& shell)
{
    WKYLIB::Mesh::readPoints(g_innerSamplePath, shell.mInnerShell);
    WKYLIB::Mesh::readPoints(g_outerSamplePath, shell.mOuterShell);
    shell.buildKdTree();
}

int main(int argc, char**argv)
{
    parseCmdLines(argc, argv);
    genDefaultParams();

    if(g_genShell) {
        SBV::Mesh mesh;
        if(jtf::mesh::load_obj(g_inputMeshPath.c_str(), mesh.triangles, mesh.vertices))
        {
            std::cout << "Fail to load mesh: " + g_inputMeshPath << std::endl;
            return 0;
        }

        auto& logger = SBV::Logger::getInstance();
        logger.setFile(g_outputPath + "/log_shell.txt");

        SBV::Shell shell;
        generateShells(mesh, shell);
    }
    else {
        SBV::Shell shell;
        readShell(shell);
        SBV::Simplifier simplifier(shell);
        simplifier.setOutputDirectory(g_outputPath);
        simplifier.setAlpha(g_alpha);
        simplifier.setGenTempResult(g_genTempResult);
        simplifier.simplify();
    }

    return 0;
}
