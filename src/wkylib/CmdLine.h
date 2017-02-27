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

        CmdLine(int argc, char** argv);
        virtual ~CmdLine();

        void parse();

        void addParamDef(const std::string& paramName, CmdParamType paramType);

        bool hasParam(const std::string& paramName) const;

        int getInt(const std::string& paramName) const;
        bool getBool(const std::string& paramName) const;
        float getFloat(const std::string& paramName) const;
        double getDouble(const std::string& paramName) const;
        const std::string& getString(const std::string& paramName) const;

    private:
        int mArgc;
        char** mArgv;

        std::map<std::string, CmdParamType> mParamDefs;
        Bundle* mBundleResult = nullptr;
    };
}

#endif
