#ifndef WKY_DEBUG_UTIL_H
#define WKY_DEBUG_UTIL_H

#include <time.h>
#include <string>
#include <exception>

namespace WKYLIB {
    class DebugTimer{
    private:
        std::string name;
        unsigned long elapsed_millisecond = 0;

        timespec time1, time2;

        bool suspended = true;

    public:
        DebugTimer(const std::string& name) : name(name) {}
        ~DebugTimer() {}

        void start()
        {
            clock_gettime(CLOCK_MONOTONIC, &time1);
            elapsed_millisecond = 0;
            suspended = false;
        }

        void suspend()
        {
            if(suspended == true)
            {
                return;
            }

            clock_gettime(CLOCK_MONOTONIC, &time2);
            elapsed_millisecond += time2.tv_sec * 1000 + time2.tv_nsec / 1000000
                    - (time1.tv_sec * 1000 + time1.tv_nsec / 1000000);
            suspended = true;
        }

        void resume()
        {
            if(suspended == false)
            {
                return;
            }
            clock_gettime(CLOCK_MONOTONIC, &time1);
            suspended = false;
        }

        void end()
        {
            if(suspended == false)
            {
                suspend();
            }
            std::cout << name << " finished after " << elapsed_millisecond << " ms." << std::endl;
        }
    };
}

#endif
