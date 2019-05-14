#ifndef __blockedarray_instance_h
#define __blockedarray_instance_h

#include <stdio.h>
#include <stdlib.h>
#include "iso3d_define.h"
#define Type Node
#define BlockedArray_Type BlockedArray_Node
#include "blockedarray.h"
#define LOG2_NODE_PER_BLOCK 12

#define Type Elem
#define BlockedArray_Type BlockedArray_Elem
#include "blockedarray.h"
#define LOG2_ELEM_PER_BLOCK 13

#endif
