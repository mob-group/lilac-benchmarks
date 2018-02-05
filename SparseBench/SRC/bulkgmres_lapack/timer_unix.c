#ifdef PentiumCPS
#include <sys/time.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>

#define CPS (PentiumCPS*1E6)
#ifndef CPS
   #define CPS (150*1E6)
#endif


static unsigned usec, sec;
static unsigned tusec, tsec;
static long long foo;

static inline void microtime(unsigned *lo, unsigned *hi)
{
  __asm __volatile (
        ".byte 0x0f; .byte 0x31   # RDTSC instruction
        movl    %%edx,%0          # High order 32 bits
        movl    %%eax,%1          # Low order 32 bits"
                : "=g" (*hi), "=g" (*lo) :: "eax", "edx");
}

double time00(void)
{
again:
  microtime(&tusec, &tsec);
  microtime(&usec, &sec);
  if (tsec != sec) goto again;

  foo = sec;
  foo = foo << 32;
  foo |= usec;
  return ((double)foo/(double)CPS);
}
#else

#include <sys/time.h>
#ifndef UseTimes
#include <sys/resource.h>
#endif

double time00(void)
{
#ifdef WALL
   struct timeval tp;
   gettimeofday(&tp, NULL);
   return( (double) (tp.tv_sec + tp.tv_usec/1000000.0) ); /* wall clock time */
#else
#ifdef UseTimes
#include <unistd.h>
   struct tms ts;
   static double ClockTick=0.0;

   if (ClockTick == 0.0) ClockTick = (double) sysconf(_SC_CLK_TCK);
   times(&ts);
   return( (double) ts.tms_utime / ClockTick );
#else
   struct rusage ruse;
   getrusage(RUSAGE_SELF, &ruse);
   return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0) );
#endif
#endif
}
#endif

double time00_(void)
{return time00();}

double TIME00(void)
{return time00();}
