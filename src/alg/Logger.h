#ifndef WKY_ALG_LOGGER_H
#define WKY_ALG_LOGGER_H

#include <string>
#include <fstream>

namespace SBV
{
    class Logger
    {
    public:
        static Logger& getInstance();

        ~Logger();

        void setFile(const std::string& path);
        void log(const std::string& text);

    private:
        Logger() = default;

    private:
        std::string mFilePath;
        std::ofstream mLogFile;
        bool mOpened = false;
    };
}

#endif
