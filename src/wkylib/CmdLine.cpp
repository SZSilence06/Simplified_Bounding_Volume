#include "CmdLine.h"
#include "Bundle.h"

namespace WKYLIB
{
    CmdLine::CmdLine(int argc, char **argv):
        mArgc(argc),
        mArgv(argv)
    {
        mBundleResult = new Bundle();
    }

    CmdLine::~CmdLine()
    {
        delete mBundleResult;
    }

    void CmdLine::addParamDef(const std::string &paramName, CmdParamType paramType)
    {
        mParamDefs.insert(std::make_pair(paramName, paramType));
    }

    bool CmdLine::hasParam(const std::string &paramName) const
    {
        return mBundleResult->hasKey(paramName);
    }

    int CmdLine::getInt(const std::string &paramName) const
    {
        return mBundleResult->getInt(paramName);
    }

    bool CmdLine::getBool(const std::string &paramName) const
    {
        return mBundleResult->getBool(paramName);
    }

    float CmdLine::getFloat(const std::string &paramName) const
    {
        return mBundleResult->getFloat(paramName);
    }

    double CmdLine::getDouble(const std::string &paramName) const
    {
        return mBundleResult->getDouble(paramName);
    }

    const std::string& CmdLine::getString(const std::string &paramName) const
    {
        return mBundleResult->getString(paramName);
    }

    bool CmdLine::parse()
    {
        return true;
    }
}
