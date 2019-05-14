/* Defines assert2(expr, message) to allow an error message to be printed. */

/* This file can be included multiple times with different settings of
   NDEBUG. */ 

#undef assert2
#undef __assert2

/* Make sure you know whether assertions should be compiled in the first
   place... */ 
#include "GR_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifdef NDEBUG
#ifdef __cplusplus
#define assert2(ig1, ig2) (static_cast<void>(0))
#else
#define assert2(ig1, ig2) ((void)(0))
#endif
#else

/* The layout of these macro definitions follows that in the GNU
   assert.h header file. */

#ifndef __GNUC__

#define assert2(expr, msg)  \
  (static_cast<void>((expr) ? 0 : __assert2 (expr, __FILE__, __LINE__, msg)))

#define __assert2(expr, file, lineno, msg)  \
  (printf ("\n%s:%u: failed assertion\nProbable cause: %s\n", \
	   file, lineno, msg), abort(), 0)

#else

#if defined(__STDC__) || defined (__cplusplus)

#define assert2(expr, msg)  \
  (static_cast<void>((expr) ? 0 : __assert2 (#expr, __FILE__, __LINE__, msg)))

#define __assert2(expr, file, line, msg)  \
  (printf ("\n%s:%u: failed assertion `%s'\nProbable cause: %s\n", \
	      file, line, expr, msg), abort(), 0)

#else /* no __STDC__ and not C++; i.e. -traditional.  */
///
extern int printf(); 

#define assert2(expr, msg)  \
  (static_cast<void>((expr) ? 0 : __assert2 (expr, __FILE__, __LINE__, msg)))

#define __assert2(expr, file, line, msg)  \
  (printf ("\n%s:%u: failed assertion `%s'\nProbable cause: %s\n", \
	      file, line, expr, msg), abort(), 0)

#endif /* no __STDC__ and not C++; i.e. -traditional.  */
#endif /* no __GNU__; i.e., /bin/cc.  */
#endif
