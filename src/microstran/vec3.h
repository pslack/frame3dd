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
MSTRANP_API double vec3_dot(vec3 A, vec3 B);
MSTRANP_API vec3 vec3_add(vec3 A, vec3 B);
MSTRANP_API vec3 vec3_cross(vec3 A, vec3 B);
MSTRANP_API vec3 vec3_scale(const vec3 A, double s);
MSTRANP_API vec3 vec3_norm(const vec3 A);
MSTRANP_API double vec3_mod(const vec3 A);

MSTRANP_API vec3 vec3_diff(const vec3 A, const vec3 B);

MSTRANP_API int vec3_print(FILE *f, const vec3 A);

MSTRANP_API vec3 vec3_rotate(vec3 A, vec3 axis, double theta);

/**
	Calculate the angle between two vectors, in radians.
*/
MSTRANP_API double vec3_angle(vec3 A, vec3 B);

/**
	Calculate the angle between two vectors, in radians. Also return
	the cross-product of the two vectors, useful with vec3_rotate.
*/
MSTRANP_API double vec3_angle_cross(vec3 A, vec3 B, vec3 *C);

MSTRANP_API char vec3_equal(vec3 A, vec3 B);
MSTRANP_API char vec3_equal_tol(vec3 A, vec3 B, double tol);

char vec3_isnan(const vec3 *A);

#define VEC3_CHECK_NAN(V) (vec3_isnan(&(V)) ? (fprintf(stderr,"%s = ",#V), vec3_print(stderr,V), fprintf(stderr,"\n"), assert(!vec3_isnan(&(V)))) : 0)

#define VEC3_PR(V) (fprintf(stderr,"%s = ",#V), vec3_print(stderr,V), fprintf(stderr,"\n"))

#define VEC3_ASSERT_EQUAL_TOL(X,Y,TOL) (vec3_equal_tol(X,Y,TOL) ? (VEC3_PR(X), VEC3_PR(Y), assert(vec3_equal_tol(X,Y,TOL))) : 0)

#ifdef __cplusplus
};
#endif

#endif

