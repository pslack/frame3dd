#define MSTRANP_BUILD
#include "vec3.h"

#include <math.h>

vec3 vec3_create(double x, double y, double z){
	vec3 A;
	A.x = x; A.y = y; A.z = z;
	return A;
}

vec3 vec3_cross(const vec3 a, const vec3 b){
	return vec3_create(
		a.y*b.z-a.z*b.y
		,a.z*b.x-a.x*b.z
		,a.x*b.y-a.y*b.x
	);
}

vec3 vec3_scale(const vec3 A, double s){
	return vec3_create(A.x * s, A.y * s, A.z * s);
}

double vec3_mod(const vec3 A){
	return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
}

vec3 vec3_norm(const vec3 A){
	return vec3_scale(A, 1./vec3_mod(A));
}

vec3 vec3_diff(const vec3 A, const vec3 B){
	return vec3_create(A.x - B.x, A.y - B.y, A.z - B.z);
}

int vec3_print(FILE *f, const vec3 A){
	return fprintf(f,"%f %f %f",A.x,A.y,A.z);
}

double vec3_dot(vec3 A, vec3 B){
	return A.x*B.x + A.y*B.y + A.z*B.z;
}

vec3 vec3_rotate(vec3 A, vec3 axis, double theta){
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

	fprintf(stderr,"x = %f, y = %f, z = %f\n",R.x, R.y, R.z);
	return vec3_scale(R, 1./SQ(ma));
#undef SQ
}

