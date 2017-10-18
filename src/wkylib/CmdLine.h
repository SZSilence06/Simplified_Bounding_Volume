#ifndef WKY_CMDLINE_H
#define WKY_CMDLINE_H

#include <map>
#include <string>

namespace WKYLIB {
    class Bundle;

    class CmdLine
    {
    public:
        enum CmdParamType
        {
            INT,
            BOOL,
            FLOAT,
            DOUBLE,
            STRING
        };

        enum ErrorCode
        {
            SUCCEEDED,
            INVALID_PARAM_NAME,
            NO_PARAM_VALUE,
            INVALID_PARAM_VALUE,
            VALUE_OVERFLOW,
            UNKNOWN_ERROR
        };

        CmdLine(int argc, char** argv);
        virtual ~CmdLine();

        bool parse();

        void addParamDef(const std::string& paramName, CmdParamType paramType);

        bool hasParam(const std::string& paramName) const;

        bool getInt(const std::string& paramName, int& value) const;
        bool getBool(const std::string& paramName, bool& value) const;
        bool getFloat(const std::string& paramName, float& value) const;
        bool getDouble(const std::string& paramName, double& value) const;
        bool getString(const std::string& paramName, std::string& value) const;

        ErrorCode getErrorCode() const;

        char* getErrorArg() const;

    private:
        int mArgc = 0;
        char** mArgv = nullptr;
        char* mErrorArg = nullptr;
        ErrorCode mErrorCode = SUCCEEDED;
        std::map<std::string, CmdParamType> mParamDefs;
        Bundle* mBundleResult = nullptr;

        bool checkParseValid(char* paramValue, char** endptr);
    };
}

#endif
