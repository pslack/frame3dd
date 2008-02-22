#ifndef MSTRAP_VEC3_H
#define MSTRAP_VEC3_H

#include "config.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

/**
	3D vector type used by the microstran parser
*/

typedef struct vec3_struct{
	double x, y, z;
} vec3;

MSTRANP_API vec3 vec3_create(double x, double y, double z);
vec3 vec3_cross(const vec3 A, const vec3 B);
MSTRANP_API vec3 vec3_scale(const vec3 A, double s);
MSTRANP_API vec3 vec3_norm(const vec3 A);
MSTRANP_API double vec3_mod(const vec3 A);

MSTRANP_API vec3 vec3_diff(const vec3 A, const vec3 B);

int vec3_print(FILE *f, const vec3 A);

MSTRANP_API vec3 vec3_rotate(vec3 A, vec3 axis, double theta);

#ifdef __cplusplus
};
#endif

#endif

