#ifndef PROFILE_H
#define PROFILE_H

#ifdef PROFILING
  #include <sys/time.h>

  #define PROFILING_INIT() struct timeval start_time, end_time, diff_time;
  #define PROFILING_START() gettimeofday(&start_time, NULL);
  #define PROFILING_END(output) gettimeofday(&end_time, NULL); timersub(&end_time, &start_time, &diff_time); printf("{P} %s: %0.3f\n", output, (double)diff_time.tv_sec*1000 + (double)diff_time.tv_usec/1000);
  #define PROFILING_END_T(output) gettimeofday(&end_time, NULL); timersub(&end_time, &start_time, &diff_time); printf("{P} [%d] %s: %0.3f\n", thrd, output, (double)diff_time.tv_sec*1000 + (double)diff_time.tv_usec/1000);
#else
  #define PROFILING_INIT()
  #define PROFILING_START()
  #define PROFILING_END(output)
  #define PROFILING_END_T(output)
#endif

#endif