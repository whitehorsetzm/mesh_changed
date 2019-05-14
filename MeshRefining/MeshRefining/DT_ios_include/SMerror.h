#if !defined(__SM_ERROR_H)
#define __SM_ERROR_H

/*
   Defines the directory where the compiled source is located; used
   in printing error messages. Each makefile has an entry
   LOCDIR     =  thedirectory
   and bmake/common includes in CFLAGS -D__SDIR__='"${LOCDIR}"'
   which is a flag passed to the compilers.
*/
#if !defined(__SDIR__)
#define __SDIR__ "unknowndirectory/"
#endif

/*
   Defines the function where the compiled source is located; used
   in printing error messages.
*/
#if !defined(__FUNC__)
#define __FUNC__ "unknownfunction"
#endif
 
/*
     These are the generic error codes. These error codes are used
     many different places in the PETSc source code.
 
*/
#define OPTMS_MEM_ERR      55  /* unable to allocate the requested memory */
#define OPTMS_ERR_FILE_OPEN       65   /* unable to open file */
#define OPTMS_ERR_FILE_READ       66   /* unable to read from file */
#define OPTMS_ERR_FILE_WRITE      67   /* unable to write to file */
#define OPTMS_ERR_FILE_UNEXPECTED 79   /* unexpected data in file */

#endif
