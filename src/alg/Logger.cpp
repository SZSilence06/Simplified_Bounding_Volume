#include "Logger.h"

namespace SBV
{
    Logger& Logger::getInstance()
    {
        static Logger logger;
        return logger;
    }

    Logger::~Logger()
    {
        if(mOpened)
        {
            mLogFile.close();
        }
    }

    void Logger::setFile(const std::string &path)
    {
        if(mFilePath == path)
            return;

        mFilePath = path;

        if(mOpened)
        {
            mLogFile.close();
        }
        mLogFile.open(path.c_str());
        mOpened = true;
    }

    void Logger::log(const std::string &text)
    {
        if(mOpened)
        {
            mLogFile << text << std::endl;
        }
    }
}
