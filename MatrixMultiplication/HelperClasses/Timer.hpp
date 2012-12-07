#ifndef TIMER_HPP
#define TIMER_HPP

#include <iostream>
#include <time.h>

class Timer {
private:
    timespec t0, t1;
    float diff;
public:
    Timer();
    Timer(bool startnow) { start(); }
    void start() { clock_gettime(CLOCK_REALTIME, &t0); }
    void stop() {
        clock_gettime(CLOCK_REALTIME, &t1);
        timespec_diff_ns(t0, t1);
    }
    float elapsed_ns() { return diff; }
    float elapsed_s() { return diff*1e-9; }
    void timespec_diff_ns(timespec &t0, timespec &t1) {
        diff = (t1.tv_sec - t0.tv_sec)*1e9 + (t1.tv_nsec - t0.tv_nsec);
    }
};

#endif
