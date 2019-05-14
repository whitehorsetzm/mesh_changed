#ifndef SUMAA_function_h
#define SUMAA_function_h

#include "SUMAA_config.h"
/* logging functions */

#if defined(__cplusplus)
   extern "C" {
#endif
void             SUMAAlogInit(int argc,char **argv);
void             SUMAAlogEventRegister(int event_num, const char *event_name);
void             SUMAAlogPrintSingle(void);
void             SUMAAlogFree(void);
double           SUMAAlogTime(void);
#ifdef PARALLEL_LOG
void             SUMAAlogPrintParallel(int nprocs, int myid, MPI_Comm procset);
#endif
#if defined(__cplusplus)
   }
#endif

#endif

