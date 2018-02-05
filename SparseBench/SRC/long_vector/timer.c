#include <stdlib.h>
#include <stdio.h>

/* Timer implementation in terms of unix rusage timer */

extern double time00(void);
double starttimer(void) {return time00();}
double starttimer_(void) {return time00();}
double STARTTIMER(void) {return time00();}
double stoptimer(void) {return time00();}
double stoptimer_(void) {return time00();}
double STOPTIMER(void) {return time00();}
