#ifndef MSTRANP_CTRANS_H
#define MSTRANP_CTRANS_H

#ifdef __cplusplus
extern "C"{
#endif

#include "config.h"
#include "vec3.h"


/**
	Coordinate transformation matrix. We can use this for translation,
	scaling, and rotation, or any combination of these. The transformation
	is stored as a 4x4 matrix. 

	This is a fairly naive implementation, so FIXME. Perhaps we can use a 3x4
	matrix instead, and perhaps we can make use of quaternions also.

	http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	http://www.win.tue.nl/~gino/solid/solid2.html
*/
typedef struct ctrans_matrix_struct{
	double m[4][4];
} ctrans_matrix;

ctrans_matrix ctrans_identity(void);

ctrans_matrix ctrans_scane(ctrans_matrix c, double s);

ctrans_matrix ctrans_rotation_z(double theta);

/**
	Apply one transformation after another, creating a new composite.
*/
ctrans_matrix ctrans_mult(ctrans_matrix A, ctrans_matrix B);

ctrans_matrix ctrans_translation(vec3 A);

/**
	Create transformation matrix for general axis/angle rotation.
*/
ctrans_matrix ctrans_rotation(vec3 axis, double theta);

/**
	Create rotation transform given the direction of the transformed Z
	and X axes.

	@return transformation that if applied to a vector in local coordinates
	relative to vectors X, Y, Z (Y = Z × X) will yeild the vector's value in
	global coordinates.
*/
MSTRANP_API ctrans_matrix ctrans_rotation_axes(vec3 Z, vec3 X);

/**
	Calculate the transform inverse (this is just a 4×4 matrix inverse)
*/
ctrans_matrix ctrans_inverse(ctrans_matrix c);

/**
	Calculate the determinant of the transform matrix.
*/
double ctrans_det(ctrans_matrix c);

/**
	Apply the coordinate transform in a 4x4 transformation matrix
	to an input vector p to return a transformed vector q
	@return the transformed vector.

	Transform:  {q} = [T]{p}.
	   | q0 |   | t00 t01 t02 t03 |   | p0 |  (x coordinate)
	   | q1 | = | t10 t11 t12 t13 | * | p1 |  (y)
	   | q2 |   | t20 t21 t22 t23 |   | p2 |  (z)
	   | 1. |   | 0.0 0.0 0.0 1.0 |   | 1. | 

	(Constants in row 4 are assumed; i.e., no perspective transformation.)
*/
MSTRANP_API vec3 ctrans_apply(ctrans_matrix c, const vec3 p);

MSTRANP_API int ctrans_print(FILE *f, const ctrans_matrix *c);

char ctrans_isnan(const ctrans_matrix *c);

#define CTRANS_CHECK_NAN(C) (ctrans_isnan(C) ? (ctrans_print(stderr,C), assert(!ctrans_isnan(C))) : 0)

#ifdef __cplusplus
};
#endif

#endif

