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
#include "ctrans.h"
#include "vec3.h"

#include <math.h>
#include <assert.h>

#define PI 3.14159265358

ctrans_matrix ctrans_identity(void){
	ctrans_matrix c;
	int i,j;
	for(i=0;i<4; ++i){
		for(j=0;j<4;++j){
			c.m[i][j] = 0.0;
		}
	};
	for(i=0;i<4;++i){
		c.m[i][i] = 1.0;
	}

	CTRANS_CHECK_NAN(&c);

	return c;
}

ctrans_matrix ctrans_scale(ctrans_matrix c, double s){
	unsigned i,j;
	for(i=0;i<4;++i){
		for(j=0;j<4;++j){
			c.m[i][j] *= s;
		}
	}
	return c;
}
	
ctrans_matrix ctrans_translation(vec3 A){
	ctrans_matrix c;
	c = ctrans_identity();
	c.m[0][3] = A.x;
	c.m[1][3] = A.y;
	c.m[2][3] = A.z;
	return c;
}

ctrans_matrix ctrans_rotation_z(double theta){
	double C = cos(theta);
	double S = sin(theta);
	ctrans_matrix c = ctrans_identity();
	c.m[0][0] = C;
	c.m[0][1] = -S;
	c.m[1][0] = S;
	c.m[1][1] = C;
	return c;
}

ctrans_matrix ctrans_mult(ctrans_matrix A, ctrans_matrix B){
	// this is just standard matrix multiplication, done naïvely
	unsigned i,j,k;
	ctrans_matrix c;
	for(i=0;i<4;++i){
		for(j=0;j<4;++j){
			double x = 0;
			for(k=0;k<4;++k){
				x += A.m[i][k] * B.m[k][j];
			}
			c.m[i][j] = x;
		}
	}
	return c;
}

char ctrans_isnan(const ctrans_matrix *c){
	/* check for NaNs */
	unsigned i,j;
	for(i=0;i<4;++i){
		for(j=0;j<4;++j){
			if(isnan(c->m[i][j])){
				return 1;
			}
		}
	}
	return 0;
}

ctrans_matrix ctrans_rotation(vec3 axis, double theta){
#if 0
	/* this is probably more efficient code, but it seems to contain a bug */
	ctrans_matrix c;

    double C = cos(theta);
    double S = sin(theta);
    double t = 1.0 - C;

	double magnitude = sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
	assert(magnitude!=0);
	if(fabs(magnitude)-1 >1e-20){
		axis.x /= magnitude;
		axis.y /= magnitude;
		axis.z /= magnitude;
	}

	C = ctrans_identity();
    C.m[0][0] = C + axis.x*axis.x*t;
    C.m[1][1] = C + axis.y*axis.y*t;
    C.m[2][2] = C + axis.z*axis.z*t;

    double tmp1 = axis.x*axis.y*t;
    double tmp2 = axis.z*S;
    C.m[1][0] = tmp1 + tmp2;
    C.m[0][1] = tmp1 - tmp2;
    tmp1 = axis.x*axis.z*t;
    tmp2 = axis.y*S;
    C.m[2][0] = tmp1 - tmp2;
    C.m[0][2] = tmp1 + tmp2;
    tmp1 = axis.y*axis.z*t;
    tmp2 = axis.x*S;
    C.m[2][1] = tmp1 + tmp2;
    C.m[1][2] = tmp1 - tmp2;
	return C;
#else
	//assert(vec3_mod(axis) > 1e-3);

	//fprintf(stderr,"C = %f, S = %f\naxis = ",C,S);
	//vec3_print(stderr,axis);
	//fprintf(stderr,"\n");
#define SQ(X) ((X)*(X))
	ctrans_matrix c;
	c = ctrans_identity();
	axis = vec3_norm(axis);

    double C = cos(theta);
    double S = sin(theta);

	assert(!isnan(C));
	assert(!isnan(S));
	assert(!isnan(axis.x));
	assert(!isnan(axis.y));
	assert(!isnan(axis.z));

	CTRANS_CHECK_NAN(&c);

//	fprintf(stderr,"c[0][1] = %f + %f = %f\n",-axis.z*S, (1-C)*axis.x*axis.y, -axis.z*S+(1-C)*axis.x*axis.y);

	// SOURCE: http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/index.htm
	c.m[0][0] = 1 + (1-C)*(SQ(axis.x)-1);
	c.m[0][1] = -axis.z*S+(1-C)*axis.x*axis.y;
	c.m[0][2] = axis.y*S+(1-C)*axis.x*axis.z;
	c.m[1][0] = axis.z*S+(1-C)*axis.x*axis.y;
	c.m[1][1] = 1 + (1-C)*(SQ(axis.y)-1);
	c.m[1][2] = -axis.x*S+(1-C)*axis.y*axis.z;
	c.m[2][0] = -axis.y*S+(1-C)*axis.x*axis.z;
	c.m[2][1] = axis.x*S+(1-C)*axis.y*axis.z;
	c.m[2][2] = 1 + (1-C)*(SQ(axis.z)-1);
#undef SQ

	CTRANS_CHECK_NAN(&c);

#if 0
	fprintf(stderr,"c[0][1] = %f = %f\n",-axis.z*S+(1-C)*axis.x*axis.y, c.m[0][1]);
	fprintf(stderr,"c[0][2] = %f = %f\n",-axis.z*S+(1-C)*axis.x*axis.y, c.m[0][1]);
	fprintf(stderr,"axis = ");
	vec3_print(stderr,axis);
	fprintf(stderr,"\ntheta = %f deg\nROT =",theta *180./PI);
	ctrans_print(stderr,&c);
#endif

#if 0
	vec3 X = vec3_create(1,0,0);

# if 0
	fprintf(stderr,"c[0][1] = %f + %f = %f\n",-axis.z*S, (1-C)*axis.x*axis.y, -axis.z*S+(1-C)*axis.x*axis.y);	ctrans_print(stderr,&c);
	fprintf(stderr,"ctrans_apply(c,X) = ");
	vec3_print(stderr,ctrans_apply(c,X));
	fprintf(stderr,"\nvec3_rotate(X,axis,theta) = ");
	vec3_print(stderr,vec3_rotate(X,axis,theta));
	fprintf(stderr,"\n");
# endif

	assert(vec3_equal_tol(ctrans_apply(c,X), vec3_rotate(X,axis,theta),1e-5));
	vec3 Y = vec3_create(0,1,0);
	assert(vec3_equal_tol(ctrans_apply(c,Y), vec3_rotate(Y,axis,theta),1e-5));
	vec3 Z = vec3_create(0,0,1);
	assert(vec3_equal_tol(ctrans_apply(c,Z), vec3_rotate(Z,axis,theta),1e-5));
	vec3 M = vec3_create(-1,-2,-3);
	assert(vec3_equal_tol(ctrans_apply(c,M), vec3_rotate(M,axis,theta),1e-5));
#endif

	return c;
}
#endif

ctrans_matrix ctrans_rotation_axes(vec3 Z, vec3 X){
	/* I was thinking about this all the wrong way! :-) */
	ctrans_matrix c = ctrans_identity();
	Z = vec3_norm(Z);
	X = vec3_norm(X);
	assert(vec3_dot(X,Z)<1e-8);
	vec3 Y = vec3_cross(Z,X);
	c.m[0][0] = X.x;
	c.m[1][0] = X.y;
	c.m[2][0] = X.z;
	c.m[0][1] = Y.x;
	c.m[1][1] = Y.y;
	c.m[2][1] = Y.z;
	c.m[0][2] = Z.x;
	c.m[1][2] = Z.y;
	c.m[2][2] = Z.z;
	return c;
}

#if 0
ctrans_matrix ctrans_rotation_axes(vec3 Z, vec3 X){
	//fprintf(stderr,"Starting ctrans_rotation_axes...\n");

	VEC3_CHECK_NAN(X);
	VEC3_CHECK_NAN(Z);

	X = vec3_norm(X);
	Z = vec3_norm(Z);

	assert(fabs(vec3_dot(X,Z))<1e-8);

	//VEC3_PR(X);
	//VEC3_PR(Z);

#define CTRANS_ROT_VEC_TOL 1e-16

	ctrans_matrix cz = ctrans_identity();
	ctrans_matrix cx = ctrans_identity();

	//fprintf(stderr,"\nRotation to orient Z...\n\n");

	vec3 z = vec3_create(0,0,1);
	vec3 xd = X;

	if(vec3_equal_tol(Z,z,CTRANS_ROT_VEC_TOL)){
		//fprintf(stderr,"No Z-rotation required\n");
	}else{
		vec3 dirz;

		assert(VEC3_NOT_EQUAL(z,Z));
		double thetaz = vec3_angle_cross(Z,z,&dirz);
		//fprintf(stderr,"Z-Rotation by %f around %f,%f,%f\n",thetaz*180./PI,dirz.x,dirz.y,dirz.z);

		/* and dirz must be normal to Z */
		assert(vec3_dot(dirz,Z)<1e-8);

		assert(vec3_mod(dirz) > 1e-3);
		cz = ctrans_rotation(dirz,thetaz);
		//CTRANS_PR(cz);

		/* apply the cz to X, this should bring X into the global x-y plane */
		xd = ctrans_apply(cz, X);
		//VEC3_PR(xd);
	}

	/* check that we recover z when applying cz to Z */
	assert(vec3_mod(vec3_diff(ctrans_apply(cz,Z),z))<1e-8);

	//fprintf(stderr,"\nRotation to orient xd to x...\n\n");
	//VEC3_PR(xd);

	assert(fabs(vec3_dot(xd,vec3_create(0,0,1)))<1e-8);

	VEC3_CHECK_NAN(xd);

	vec3 dirx;
	VEC3_CHECK_NAN(dirx);

	vec3 x = vec3_create(1,0,0);
	if(vec3_equal_tol(xd,x,CTRANS_ROT_VEC_TOL)){
		//fprintf(stderr,"No X rotation required\n");
	}else{
		double thetax = 0;
		//VEC3_PR(xd);
		//VEC3_PR(x);
		thetax = vec3_angle_cross(xd,x,&dirx);
		//fprintf(stderr,"thetax = %f\n",thetax*180./PI);
		VEC3_CHECK_NAN(dirx);
		if(vec3_mod(dirx)<1e-8){
			//fprintf(stderr,"dirx is very small!\n");
			/* this means that xd × x is close to zero, hence theta is close 
			to 0 or ±180. so we'll use 0,0,1 as the axis in that case, becase
			the orientation 0,0,±1 won't matter */
			dirx = vec3_create(0,0,1);
			
		}
		//fprintf(stderr,"X-Rotation by %f around %f,%f,%f\n",thetax*180./PI,dirx.x,dirx.y,dirx.z);
		//assert(vec3_mod(dirx)>1e-3);	
		cx = ctrans_rotation(dirx,thetax);
	}
#if 1

	VEC3_ASSERT_EQUAL_TOL(vec3_create(0,0,1),vec3_create(0,0,1),1e-8);
	//VEC3_ASSERT_EQUAL_TOL(vec3_create(0,0,1),vec3_create(0,0,1.001),1e-8);

	ctrans_matrix c = ctrans_mult(cx,cz);
	//CTRANS_PR(c);

	// test the transform: 'X' must transfor to 'x'...
	vec3 gX = ctrans_apply(c,X);
	VEC3_ASSERT_EQUAL_TOL(gX,vec3_create(1,0,0),1e-8);

	// and 'x' must transform to 'X'...
	ctrans_matrix ci = ctrans_inverse(c);
	gX = ctrans_apply(ci,x);
	VEC3_ASSERT_EQUAL_TOL(gX,X,1e-8);
	
#endif
	
	return ctrans_mult(cx, cz);
}
#endif

vec3 ctrans_apply(ctrans_matrix c, const vec3 p){
	vec3 q;
	q.x = c.m[0][0] * p.x + c.m[0][1] * p.y + c.m[0][2] * p.z + c.m[0][3];
	q.y = c.m[1][0] * p.x + c.m[1][1] * p.y + c.m[1][2] * p.z + c.m[1][3];
	q.z = c.m[2][0] * p.x + c.m[2][1] * p.y + c.m[2][2] * p.z + c.m[2][3];
	return q;
}

int ctrans_print(FILE *f, const ctrans_matrix *c){
	int n = 0;
	unsigned i;
	for(i=0; i<4;++i){
		n += fprintf(f,"\t[ %8f %8f %8f %8f ]\n",c->m[i][0],c->m[i][1],c->m[i][2],c->m[i][3]);
	}
	return n;
}

MSTRANP_API char ctrans_equal_tol(const ctrans_matrix *c, const ctrans_matrix *d, double tol){
	unsigned i,j;
	for(i=0;i<4;++i){
		for(j=0;j<4;++j){
			if(fabs(c->m[i][j] - d->m[i][j]) > tol)return 0;
		}
	}
	return 1;
}


ctrans_matrix ctrans_inverse(ctrans_matrix c){
#define m00 c.m[0][0]
#define m01 c.m[0][1]
#define m02 c.m[0][2]
#define m03 c.m[0][3]
#define m10 c.m[1][0]
#define m11 c.m[1][1]
#define m12 c.m[1][2]
#define m13 c.m[1][3]
#define m20 c.m[2][0]
#define m21 c.m[2][1]
#define m22 c.m[2][2]
#define m23 c.m[2][3]
#define m30 c.m[3][0]
#define m31 c.m[3][1]
#define m32 c.m[3][2]
#define m33 c.m[3][3]
	ctrans_matrix v;
	v.m[0][0] = m12*m23*m31 - m13*m22*m31 + m13*m21*m32 - m11*m23*m32 - m12*m21*m33 + m11*m22*m33;
	v.m[0][1] = m03*m22*m31 - m02*m23*m31 - m03*m21*m32 + m01*m23*m32 + m02*m21*m33 - m01*m22*m33;
	v.m[0][2] = m02*m13*m31 - m03*m12*m31 + m03*m11*m32 - m01*m13*m32 - m02*m11*m33 + m01*m12*m33;
	v.m[0][3] = m03*m12*m21 - m02*m13*m21 - m03*m11*m22 + m01*m13*m22 + m02*m11*m23 - m01*m12*m23;
	v.m[1][0] = m13*m22*m30 - m12*m23*m30 - m13*m20*m32 + m10*m23*m32 + m12*m20*m33 - m10*m22*m33;
	v.m[1][1] = m02*m23*m30 - m03*m22*m30 + m03*m20*m32 - m00*m23*m32 - m02*m20*m33 + m00*m22*m33;
	v.m[1][2] = m03*m12*m30 - m02*m13*m30 - m03*m10*m32 + m00*m13*m32 + m02*m10*m33 - m00*m12*m33;
	v.m[1][3] = m02*m13*m20 - m03*m12*m20 + m03*m10*m22 - m00*m13*m22 - m02*m10*m23 + m00*m12*m23;
	v.m[2][0] = m11*m23*m30 - m13*m21*m30 + m13*m20*m31 - m10*m23*m31 - m11*m20*m33 + m10*m21*m33;
	v.m[2][1] = m03*m21*m30 - m01*m23*m30 - m03*m20*m31 + m00*m23*m31 + m01*m20*m33 - m00*m21*m33;
	v.m[2][2] = m01*m13*m30 - m03*m11*m30 + m03*m10*m31 - m00*m13*m31 - m01*m10*m33 + m00*m11*m33;
	v.m[2][3] = m03*m11*m20 - m01*m13*m20 - m03*m10*m21 + m00*m13*m21 + m01*m10*m23 - m00*m11*m23;
	v.m[3][0] = m12*m21*m30 - m11*m22*m30 - m12*m20*m31 + m10*m22*m31 + m11*m20*m32 - m10*m21*m32;
	v.m[3][1] = m01*m22*m30 - m02*m21*m30 + m02*m20*m31 - m00*m22*m31 - m01*m20*m32 + m00*m21*m32;
	v.m[3][2]= m02*m11*m30 - m01*m12*m30 - m02*m10*m31 + m00*m12*m31 + m01*m10*m32 - m00*m11*m32;
	v.m[3][3] = m01*m12*m20 - m02*m11*m20 + m02*m10*m21 - m00*m12*m21 - m01*m10*m22 + m00*m11*m22;
	return ctrans_scale(v, 1./ctrans_det(c));
}

double ctrans_det(ctrans_matrix c){
	return 
	   m03 * m12 * m21 * m30-m02 * m13 * m21 * m30-m03 * m11 * m22 * m30+m01 * m13 * m22 * m30
     + m02 * m11 * m23 * m30-m01 * m12 * m23 * m30-m03 * m12 * m20 * m31+m02 * m13 * m20 * m31
     + m03 * m10 * m22 * m31-m00 * m13 * m22 * m31-m02 * m10 * m23 * m31+m00 * m12 * m23 * m31
     + m03 * m11 * m20 * m32-m01 * m13 * m20 * m32-m03 * m10 * m21 * m32+m00 * m13 * m21 * m32
     + m01 * m10 * m23 * m32-m00 * m11 * m23 * m32-m02 * m11 * m20 * m33+m01 * m12 * m20 * m33
     + m02 * m10 * m21 * m33-m00 * m12 * m21 * m33-m01 * m10 * m22 * m33+m00 * m11 * m22 * m33;
} 

