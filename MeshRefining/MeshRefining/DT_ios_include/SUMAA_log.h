#include "SUMAA_config.h"
/*@ SUMAA_log.h - This file defines the type of logging done in the
                  SUMAA3d package.  

	System Description:
    The facility is initialized by calling the macro SUMAA_LOG_INIT.
    The facility is terminated and statistics printed by calling the
    macro SUMAA_LOG_PRINT.  An new event is registered by calling the macro 
    SUMAA_LOG_EVENT_REGISTER and the logging of the event is begun by 
    calling SUMAA_LOG_EVENT_BEGIN.  At then end of an event call the
    macro SUMAA_LOG_FLOPS to record the number of flops associated with the
    event and SUMAA_LOG_EVENT_END.

    Somewhere in the users main program the preprocessor variable
    SUMAA_MAIN_LOG must be defined (#define SUMAA_MAIN_LOG) if logging is turned
    on.  This allows the declaration of certain global variables.

    To print the logging informaiton, the default is to assume a single
    processor.  For parallel logging to be printed properly, you must
    define #PARALLEL_LOG and send the procinfo structure to the logging
    print routine.

@*/
#ifndef SUMAA_LOG_INCLUDE_DEFINED
#define SUMAA_LOG_INCLUDE_DEFINED

#ifdef __cplusplus
extern "C" {
#endif

#ifdef mpi
#include "SUMAA_mpi.h"
#endif

#ifdef mpi
#define SUMAA_LOG_exit(a) MPI_abort(MPI_COMM_WORLD,a)
#else 
#define SUMAA_LOG_exit(a) exit(a)
#endif

#include <stdio.h>

#if defined(SUMAA_LOG)

#if defined(SUMAA_LOG_LITE)
/* don't do max min or nested flop counts*/
#define MAX_MIN_LOG(a) 
#define NESTED_LOG(a) 
#else
#define MAX_MIN_LOG(a) a
#define NESTED_LOG(a) a
#endif

/* define the clock to use */
#if defined (MPI_CLOCK)
#define LOG_clock() MPI_Wtime();
#else
#define LOG_clock() SUMAAlogTime();
#endif

#define SUMAA_NO_PARENT       -1
#define SUMAA_LOG_MAX_EVENTS   300
#define SUMAA_USER_EVENT_LOW   350
#define SUMAA_USER_EVENT_HIGH  SUMAA_USER_EVENT_LOW+50
    
/* define a structure that contains all information regarding an event */
typedef struct {
	    char    *name;
	    int     num_calls;
        int     parent;
        int     level;
	    double	time_stamp;
	    double	total_time;
        double  event_time;
	    double	max_time;
	    double	min_time;
	    double  total_flops;
	    double  event_flops;
	    double  min_flops;
	    double  max_flops;
        int     registered;
} SUMAA_LOG_struct;
    
typedef struct {
        double item;
        double sum;
        double min;
        double max;
        double avg;
} SUMAA_LOG_maxmind;

    
/* a set of global variables that are used in logging */
#ifdef SUMAA_MAIN_LOG
        int              SUMAA_LOG_current;
        int              SUMAA_LOG_initialized = 0;
        double           SUMAA_LOG_total_flops;
        double           SUMAA_LOG_time;
        char            *SUMAA_LOG_file;
        SUMAA_LOG_struct SUMAA_LOG_event[SUMAA_LOG_MAX_EVENTS];
#else
        extern int              SUMAA_LOG_current;
        extern int              SUMAA_LOG_initialized;
        extern double           SUMAA_LOG_total_flops;
        extern double           SUMAA_LOG_time;
        extern char            *SUMAA_LOG_file;
        extern SUMAA_LOG_struct SUMAA_LOG_event[SUMAA_LOG_MAX_EVENTS];
#endif
    
#else 
#define SUMAA_USER_EVENT_LOW
#endif

#ifdef __cplusplus
}
#endif

#include "SUMAA_macro.h"
#include "SUMAA_function.h"

#endif
