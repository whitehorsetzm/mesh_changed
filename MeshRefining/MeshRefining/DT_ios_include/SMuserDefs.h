#ifndef SM_USER_DEFS_H
#define SM_USER_DEFS_H

#define OPTMS_DEFAULT -1

/* 2D function/gradient options included with Opt-MS
   The default is MAX_MIN_SINE */
#define MAX_MIN_ANGLE 1
#define MIN_MAX_COSINE 2
#define MAX_MIN_COSINE 3
#define MAX_MIN_SINE 4
#define MIN_MAX_ANGLE 5
#define MIN_MAX_JACOBIAN_DIFF 6
#define MAX_MIN_SCALED_JACOBIAN 7
#define MAX_MIN_AREA_LENGTH_RATIO 8
#define MIN_MAX_LENGTH_AREA_RATIO 9
#define MAX_MIN_INTERIOR_ANGLE 10
#define MAX_MIN_INTERIOR_SINE 11
#define MIN_MAX_INTERIOR_COSINE 12
#define MAX_MIN_INTERIOR_SCALED_JACOBIAN 13
#define MIN_MAX_NORM_JAC_SQUARED_2D 14
#define FUNCTION_DEFAULT_2D 4

/* 3D function/gradient options included with Opt-MS
   The default is MAX_SINE_DIHEDRAL */
#define MAX_MIN_DIHEDRAL 21
#define MIN_MAX_DIHEDRAL 22
#define MAX_MIN_COSINE_DIHEDRAL 23
#define MIN_MAX_COSINE_DIHEDRAL 24
#define MAX_MIN_SINE_DIHEDRAL 25
#define MAX_MIN_SCALED_JACOBIAN_3D 26
#define MIN_MAX_SRMS_VOLUME_RATIO 27
#define MIN_MAX_CONDITION 28
#define MIN_MAX_NORM_JAC_SQUARED_3D 29
#define MAX_MIN_SOLID 30
#define FUNCTION_DEFAULT_3D 25

/* user interface to the smoothing techniques */
#define OPTMS_LAPLACIAN          'S'
#define OPTMS_SMART_LAPLACIAN    'L'
#define OPTMS_OPTIMIZATION_ONLY  'O'
#define OPTMS_COMBINED           'C'
#define OPTMS_COMBINED1          '1'
#define OPTMS_COMBINED2          '2'
#define OPTMS_COMBINED3          '3'
#define OPTMS_FLOATING_THRESHOLD 'F'
#define OPTMS_STUPID_LAPLACE     'S'
#define OPTMS_TECHNIQUE_DEFAULT  'C'

/* smoothing techniques */
#define LAPLACE_ONLY      0
#define OPTIMIZATION_ONLY 1
#define COMBINED          2
#define COMBINED1         3
#define COMBINED2         4
#define COMBINED3         5
#define FLOATING_THRESHOLD 6
#define STUPID_LAPLACE    7
#define DEFAULT_LAP_ACCEPT -1
#define TECHNIQUE_DEFAULT 2

/* untangling techniques */
#define LAPLACE_ONLY        0
#define LINEAR_PROGRAM_ONLY 1
#define COMBINED_UNTANGLING 2


#endif





