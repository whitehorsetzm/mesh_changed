#ifndef SUMAA_MACRO_H
#define SUMAA_MACRO_H

#include "SUMAA_config.h"
#include <assert.h>
#ifdef SUMAA_LOG
#define SUMAA_LOG_INIT(argc,argv)        SUMAAlogInit(argc,argv)
#define SUMAA_LOG_EVENT_REGISTER(a,str)  SUMAAlogEventRegister(a,str)

#ifdef PARALLEL_LOG
#define SUMAA_LOG_PRINT(a,b,c)           SUMAAlogPrintParallel(a,b,c)
#else
#define SUMAA_LOG_PRINT                  SUMAAlogPrintSingle()
#endif

#define SUMAA_LOG_FREE                   SUMAAlogFree()

#define SUMAA_LOG_EVENT_BEGIN(a)         SUMAAlogEventBegin(a)
#define SUMAA_LOG_EVENT_END(a)           SUMAAlogEventEnd(a)

#define SUMAAlogEventBegin(event_num) \
{ \
    assert(event_num < SUMAA_LOG_MAX_EVENTS); \
    if (SUMAA_LOG_event[event_num].registered == 0) { \
     /*fprintf(stderr,"Warning: Trying to log an unregistered event number %d\n",\
                     event_num);*/\
    } else { \
     SUMAA_LOG_event[event_num].num_calls++; \
     SUMAA_LOG_event[event_num].event_flops = 0; \
     SUMAA_LOG_event[event_num].event_time = 0; \
     if (SUMAA_LOG_current != event_num) { \
       SUMAA_LOG_event[event_num].parent = SUMAA_LOG_current; \
     } else { \
       /*fprintf(stderr, "Nesting:  Event #%d in %s, line %d\n", event_num, __FILE__, __LINE__);*/ \
       SUMAA_LOG_event[event_num].level++; \
     } \
     SUMAA_LOG_current = event_num; \
     SUMAA_LOG_event[event_num].time_stamp = LOG_clock(); \
    } \
}

#define SUMAA_LOG_FLOPS(a,b) \
{ \
    if (SUMAA_LOG_event[a].registered == 0) { \
     fprintf(stderr,"Warning: Trying to log unregistered event (No. %d) flops\n", a);\
    } \
    SUMAA_LOG_event[a].event_flops += b; \
}

#define SUMAA_LOG_GLOBAL_TIME(a)  \
{ \
    if (SUMAA_LOG_event[a].registered == 0) { \
     fprintf(stderr,"Warning: Trying to add unregistered event (No. %d) time to global time \n", a);\
    } \
    SUMAA_LOG_time += SUMAA_LOG_event[a].event_time; \
}

#define SUMAAlogEventEnd(event_num) \
{ \
  double time_end; \
  double event_time, event_flops; \
  assert(event_num < SUMAA_LOG_MAX_EVENTS); \
  if (SUMAA_LOG_event[event_num].registered == 0) { \
   /* fprintf(stderr,"Warning: Trying to log an unregistered event number %d\n",\
	    event_num);*/\
  } else { \
    /* Record timing information */ \
    time_end = LOG_clock(); \
    event_time = time_end - SUMAA_LOG_event[event_num].time_stamp; \
    SUMAA_LOG_event[event_num].event_time = event_time; \
    SUMAA_LOG_event[event_num].total_time += event_time; \
     \
    /* Record the flop information */ \
    event_flops = SUMAA_LOG_event[event_num].event_flops; \
    SUMAA_LOG_event[event_num].total_flops += event_flops; \
    SUMAA_LOG_total_flops += event_flops; \
    \
    MAX_MIN_LOG({ \
      /* max min time per call */ \
      if (event_time < SUMAA_LOG_event[event_num].min_time) { \
	SUMAA_LOG_event[event_num].min_time = event_time; \
      } \
      if (event_time > SUMAA_LOG_event[event_num].max_time) { \
	SUMAA_LOG_event[event_num].max_time = event_time; \
      } \
     \
      /* max min flops per call */ \
      if (event_flops < SUMAA_LOG_event[event_num].min_flops) { \
	SUMAA_LOG_event[event_num].min_flops = event_flops; \
      } \
      if (event_flops > SUMAA_LOG_event[event_num].max_flops) { \
	SUMAA_LOG_event[event_num].max_flops = event_flops; \
      } \
    }); \
     \
    NESTED_LOG({ \
      int parent; \
      int level; \
      if (SUMAA_LOG_initialized) { \
	level = SUMAA_LOG_event[event_num].level; \
	parent = SUMAA_LOG_event[event_num].parent; \
	while (parent != SUMAA_NO_PARENT) { \
	  if (level == 0) { \
	    SUMAA_LOG_event[parent].total_flops += event_flops; \
	    parent = SUMAA_LOG_event[parent].parent; \
	  } else { \
	    level--; \
	  } \
	} \
      } else { \
	fprintf(stderr,"You haven't initialized the logging routines\n"); \
	fprintf(stderr,"Not doing nested logging\n"); \
      } \
    }); \
    /* update the current event and level information */ \
    assert(SUMAA_LOG_current != SUMAA_LOG_event[event_num].parent); \
    if ((SUMAA_LOG_current != SUMAA_LOG_event[event_num].parent) &&  \
	(SUMAA_LOG_event[event_num].level == 0)) { \
      SUMAA_LOG_current = SUMAA_LOG_event[event_num].parent; \
    } else if (SUMAA_LOG_current != SUMAA_LOG_event[event_num].parent){\
      SUMAA_LOG_event[event_num].level--; \
    } \
  } \
} 

/* define the communication routines for sum/min/max 
   using mpi for the parallel logging */
#if defined(PARALLEL_LOG) && defined(mpi)

#define LOG_GDSUM(sum_vec,vec_len,work_vec,procset) \
{ \
    int i99; \
    double  *dptr1 = (double *) sum_vec, *dptr2 = (double *) work_vec; \
    for (i99=0;i99<(vec_len);i99++) { \
        (dptr2)[i99] = (dptr1)[i99]; \
    } \
    MPI_Allreduce(dptr2,dptr1,vec_len,MPI_DOUBLE,MPI_SUM,(MPI_Comm)procset); \
}

#define LOG_GDMIN(sum_vec,vec_len,work_vec,procset) \
{ \
    int i99; \
    double  *dptr1 = (double *) sum_vec, *dptr2 = (double *) work_vec; \
    for (i99=0;i99<(vec_len);i99++) { \
        (dptr2)[i99] = (dptr1)[i99]; \
    } \
    MPI_Allreduce(dptr2,dptr1,vec_len,MPI_DOUBLE,MPI_MIN,(MPI_Comm)procset); \
}

#define LOG_GDMAX(sum_vec,vec_len,work_vec,procset) \
{ \
    int i99; \
    double  *dptr1 = (double *) sum_vec, *dptr2 = (double *) work_vec; \
    for (i99=0;i99<(vec_len);i99++) { \
        (dptr2)[i99] = (dptr1)[i99]; \
    } \
    MPI_Allreduce(dptr2,dptr1,vec_len,MPI_DOUBLE,MPI_MAX,(MPI_Comm)procset); \
}
#endif

#define SUMAA_INIT_MAX_MIN(a,value) \
{ \
    a.item = (double) value; \
    a.sum = 0.;  a.max = 0.;  a.min=0.;  a.avg=0.; \
}

#define SUMAA_MIN_MAX_AVGd(a,nprocs,procset) \
{ \
    double t, sum99, max99, min99, avg99; \
    sum99 = a.item; max99 = a.item;  min99 = a.item; \
    LOG_GDSUM(&sum99,1,&t,procset); \
    LOG_GDMAX(&max99,1,&t,procset); \
    LOG_GDMIN(&min99,1,&t,procset); \
    avg99 = sum99/nprocs; \
    a.sum = sum99; a.max = max99; a.min = min99; a.avg = avg99; \
}

#else

#define SUMAA_LOG_INIT(a,b)
#define SUMAA_LOG_EVENT_REGISTER(a,str) 
#define SUMAA_LOG_FLOPS(a,b) 
#define SUMAA_LOG_GLOBAL_TIME(a)  
#define SUMAA_LOG_EVENT_BEGIN(a)
#define SUMAA_LOG_EVENT_END(a)  
#define SUMAA_LOG_PRINT       
#define SUMAA_LOG_FREE

#endif

/*  #if defined(rs6000) */
/*  #define SUMAA_LOG_time_func(v) {static struct  timestruc_t _tp; \ */
/*                                 UTP_readTime(&_tp); \ */
/*                                 (v)=((double)_tp.tv_sec)+(1.0e-9)*(_tp.tv_nsec);} */

/*  #else */
#define SUMAA_LOG_time_func(v) {static struct timeval _tp; \
                               gettimeofday(&_tp,(struct timezone *)0);\
                               (v)=((double)_tp.tv_sec)+(1.0e-6)*(_tp.tv_usec);}
/*  #endif */

#endif
