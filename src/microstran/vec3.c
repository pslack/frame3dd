/*	FRAME3DD: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 2007-2008 John Pye

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define MSTRANP_BUILD
#include "vec3.h"

#include <math.h>
#include <assert.h>

const vec3 VEC3_ZERO = {0,0,0};

vec3 vec3_create(double x, double y, double z){
	vec3 A;
	A.x = x; A.y = y; A.z = z;
	return A;
}

vec3 vec3_add(vec3 A, vec3 B){
	return vec3_create(
		A.x+B.x, A.y + B.y, A.z + B.z
	);
}

double vec3_dot(vec3 a, vec3 b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

vec3 vec3_cross(vec3 a, vec3 b){
	return vec3_create(
		a.y*b.z-a.z*b.y
		,a.z*b.x-a.x*b.z
		,a.x*b.y-a.y*b.x
	);
}

vec3 vec3_scale(vec3 A, double s){
	return vec3_create(A.x * s, A.y * s, A.z * s);
}

double vec3_mod(vec3 A){
	return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
}

vec3 vec3_norm(vec3 A){
	return vec3_scale(A, 1./vec3_mod(A));
}

vec3 vec3_diff(vec3 A, vec3 B){
	return vec3_create(A.x - B.x, A.y - B.y, A.z - B.z);
}

vec3 vec3_negate(vec3 A){
	vec3 B;
	B.x = -A.x;
	B.y = -A.y;
	B.z = -A.z;
	return B;
}

int vec3_print(FILE *f, vec3 A){
	return fprintf(f,"%.10f %.10f %.10f",A.x,A.y,A.z);
}

vec3 vec3_rotate(vec3 A, vec3 axis, double theta){
#if 0
	assert(!isnan(theta));
	assert(!vec3_isnan(&axis));
	assert(!vec3_isnan(&A));
#endif

	/* source: http://www.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html */
	double Aa = vec3_dot(A,axis);
	double ma = vec3_mod(axis);

	double u = axis.x, v = axis.y, w = axis.z;
	double ct = cos(theta);
	double st = sin(theta);
	vec3 R;
#define SQ(X) ((X)*(X))
	R.x = u*Aa + ( A.x*(SQ(v)+SQ(w)) - u*(v*A.y+w*A.z) )*ct + ma * (-w*A.y+v*A.z) * st;
	R.y = v*Aa + ( A.y*(SQ(u)+SQ(w)) - v*(u*A.x+w*A.z) )*ct + ma * (w*A.x-u*A.z) * st;
	R.z = w*Aa + ( A.z*(SQ(u)+SQ(v)) - w*(u*A.x+v*A.y) )*ct + ma * (-v*A.x+u*A.y) * st;

#if 0
	if(isnan(R.x) || isnan(R.y) || isnan(R.z)){
		VEC3_PR(R);
	}
	assert(!isnan(R.x));
	assert(!isnan(R.y));
	assert(!isnan(R.z));
#endif

	/* fprintf(stderr,"x = %f, y = %f, z = %f\n",R.x, R.y, R.z); */
	return vec3_scale(R, 1./SQ(ma));
#undef SQ
}

double vec3_angle(vec3 A, vec3 B){
	assert(A.x!=B.x || A.y!=B.y || A.z!=B.z);
	vec3 C = vec3_cross(A,B);
	return atan2(vec3_mod(C),vec3_dot(A,B));
}

MSTRANP_API double vec3_angle_cross(vec3 A, vec3 B, vec3 *C){
	assert(VEC3_NOT_EQUAL(A,B));
	VEC3_CHECK_NAN(A);
	VEC3_CHECK_NAN(B);
	assert(C!=NULL);
	C->x = A.y*B.z-A.z*B.y;
	C->y = A.z*B.x-A.x*B.z;
	C->z = A.x*B.y-A.y*B.x;
	if(vec3_isnan(C)){
		VEC3_PR(A);
		VEC3_PR(B);
		//fprintf(stderr,"dot(A,B) = %f\n",vec3_dot(A,B));
		VEC3_CHECK_NAN(*C);
	}
	return atan2(vec3_mod(*C),vec3_dot(A,B));
}

char vec3_equal_tol(vec3 A, vec3 B, double tol){
	return vec3_mod(vec3_diff(B,A))<tol;
}

char vec3_equal(vec3 A, vec3 B){
	return A.x==B.x && A.y==B.y && A.z==B.z;
}

char vec3_isnan(const vec3 *A){
	return isnan(A->x) || isnan(A->y) || isnan(A->z);
}


