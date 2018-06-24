#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef __MACH__
   #include <mach/clock.h>
   #include <mach/mach.h>
#endif

double fsils_hrcput_()
{
   struct timespec time;
   double res;

   #ifdef __MACH__    // OSX does not have clock_gettime
      clock_serv_t cclock;
      mach_timespec_t mts;
      host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
      clock_get_time(cclock, &mts);
      mach_port_deallocate(mach_task_self(), cclock);
      time.tv_sec = mts.tv_sec;
      time.tv_nsec = mts.tv_nsec;
   #else
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time);
   #endif

   res = ((double)time.tv_nsec)/1000000000.0 + (double)(time.tv_sec%1000000);

   return res;
}
 
