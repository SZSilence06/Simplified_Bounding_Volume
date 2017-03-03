#include "CmdLine.h"
#include "Bundle.h"
#include <string.h>
#include <cerrno>

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

    bool CmdLine::getInt(const std::string &paramName, int& value) const
    {
        return mBundleResult->getInt(paramName, value);
    }

    bool CmdLine::getBool(const std::string &paramName, bool& value) const
    {
        return mBundleResult->getBool(paramName, value);
    }

    bool CmdLine::getFloat(const std::string &paramName, float& value) const
    {
        return mBundleResult->getFloat(paramName, value);
    }

    bool CmdLine::getDouble(const std::string &paramName, double& value) const
    {
        return mBundleResult->getDouble(paramName, value);
    }

    bool CmdLine::getString(const std::string &paramName, std::string& value) const
    {
        return mBundleResult->getString(paramName, value);
    }

    bool CmdLine::parse()
    {
        char* paramName = nullptr;
        char* paramValue = nullptr;

        for(int i = 1; i < mArgc; i++)
        {
            paramName = mArgv[i];

            auto iter = mParamDefs.find(paramName);
            if(iter == mParamDefs.end())
            {
                mErrorCode = INVALID_PARAM_NAME;
                mErrorArg = paramName;
                return false;
            }
            CmdParamType paramType = iter->second;
            if(paramType == CmdParamType::BOOL)
            {
                mBundleResult->putBool(paramName, true);
                continue;
            }
            if(++i >= mArgc)
            {
                mErrorCode = NO_PARAM_VALUE;
                mErrorArg = paramName;
                return false;
            }

            paramValue = mArgv[i];
            switch(paramType)
            {
            case CmdParamType::INT:
            {
                char **endptr;
                int value = std::strtol(paramValue, endptr, 10);
                if(checkParseValid(paramValue, endptr) == false)
                {
                    return false;
                }
                mBundleResult->putInt(paramName, value);
                break;
            }
            case CmdParamType::FLOAT:
            {
                char **endptr;
                float value = std::strtof(paramValue, endptr);
                if(checkParseValid(paramValue, endptr) == false)
                {
                    return false;
                }
                mBundleResult->putFloat(paramName, value);
                break;
            }
            case CmdParamType::DOUBLE:
            {
                char **endptr;
                double value = std::strtod(paramValue, endptr);
                if(checkParseValid(paramValue, endptr) == false)
                {
                    return false;
                }
                mBundleResult->putDouble(paramName, value);
                break;
            }
            case CmdParamType::STRING:
            {
                mBundleResult->putString(paramName, paramValue);
                break;
            }
            default:
                mErrorCode = UNKNOWN_ERROR;
                return false;
            }
        }

        mErrorCode = SUCCEEDED;
        return true;
    }

    CmdLine::ErrorCode CmdLine::getErrorCode() const
    {
        return mErrorCode;
    }
    char* CmdLine::getErrorArg() const
    {
        return mErrorArg;
    }

    bool CmdLine::checkParseValid(char *paramValue, char **endptr)
    {
        if(*endptr == paramValue)
        {
            mErrorCode = INVALID_PARAM_VALUE;
            mErrorArg = paramValue;
            return false;
        }
        if(errno == ERANGE)
        {
            mErrorCode = VALUE_OVERFLOW;
            mErrorArg = paramValue;
            return false;
        }
        return true;
    }
}
