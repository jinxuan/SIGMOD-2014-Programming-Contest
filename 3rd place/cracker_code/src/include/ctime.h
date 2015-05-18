/*
 * Author: Shao Yingxia
 * Create Date: 2012年12月12日 星期三 22时55分34秒
 */
#ifndef __CTIME_H__
#define __CTIME_H__

#include <sys/time.h>

class Time {
    public:
        Time();
        ~Time();
        void Start();
        void Reset();
        void Stop();
        double GetElapsedTime();
    private:
        struct timeval start_;
        struct timeval end_;
        struct timezone tz_;
        bool isstart_;
};

#endif    // #ifndef __TIME_H__

