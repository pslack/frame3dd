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
	double C = cos(theta);
	double S = sin(theta);
	fprintf(stderr,"C = %f, S = %f\naxis = ",C,S);
	vec3_print(stderr,axis);
	fprintf(stderr,"\n");
#define SQ(X) ((X)*(X))
	ctrans_matrix c;
	c = ctrans_identity();
	axis = vec3_norm(axis);

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

ctrans_matrix ctrans_rotation_axes(vec3 Z, vec3 X){

	VEC3_CHECK_NAN(X);
	VEC3_CHECK_NAN(Z);

	X = vec3_norm(X);
	Z = vec3_norm(Z);

	VEC3_PR(X);
	VEC3_PR(Z);

#define CTRANS_ROT_VEC_TOL 1e-8

	ctrans_matrix cz = ctrans_identity();
	ctrans_matrix cx = ctrans_identity();

	fprintf(stderr,"\nRotation to orient Z...\n\n");

	vec3 z = vec3_create(0,0,1);
	vec3 xd = X;

	if(!vec3_equal_tol(Z,z,CTRANS_ROT_VEC_TOL)){
		fprintf(stderr,"No Z-rotation required\n");
	}else{
		vec3 dirz;
		double thetaz = vec3_angle_cross(z,Z,&dirz);
		fprintf(stderr,"Z-Rotation by %f around %f,%f,%f\n",thetaz*180./PI,dirz.x,dirz.y,dirz.z);

		/* and dirz must be normal to Z */
		assert(vec3_dot(dirz,Z)<1e-8);

		cz = ctrans_rotation(dirz,-thetaz);
		CTRANS_PR(cz);

		/* apply the cz to X, this should bring X into the global x-y plane */
		vec3 xd = ctrans_apply(cz, X);
	}

	/* check that we recover z when applying cz to Z */
	assert(vec3_mod(vec3_diff(ctrans_apply(cz,Z),z))<1e-8);

	fprintf(stderr,"\nRotation to orient X...\n\n");

	assert(fabs(vec3_dot(xd,vec3_create(0,0,1)))<1e-8);

	VEC3_CHECK_NAN(xd);

	VEC3_PR(xd);
	vec3 dirx;
	VEC3_CHECK_NAN(dirx);

	vec3 x = vec3_create(1,0,0);
	if(vec3_equal_tol(xd,x,1e-8)){
		fprintf(stderr,"No X rotation required\n");
	}else{
		double thetax = 0;
		thetax = vec3_angle_cross(xd,x,&dirx);
		fprintf(stderr,"thetax = %f\n",thetax);
		VEC3_CHECK_NAN(dirx);
		if(thetax<1e-5){
			fprintf(stderr,"No X rotation required\n");
			return cz;
		}
		fprintf(stderr,"X-Rotation by %f around %f,%f,%f\n",thetax*180./PI,dirx.x,dirx.y,dirx.z);
		cx = ctrans_rotation(dirx,-thetax);
	}
#if 1

	VEC3_ASSERT_EQUAL_TOL(vec3_create(0,0,1),vec3_create(0,0,1),1e-8);
	VEC3_ASSERT_EQUAL_TOL(vec3_create(0,0,1),vec3_create(0,0,1.001),1e-8);

	// test the transform
	ctrans_matrix c = ctrans_mult(cx,cz);
	ctrans_matrix ci = ctrans_inverse(c);
	vec3 gX = ctrans_apply(c,vec3_create(1,0,0));
	VEC3_ASSERT_EQUAL_TOL(gX,X,1e-8);

	gX = ctrans_apply(ci,gX);
	vec3 xxx = vec3_create(1,0,0);
	VEC3_ASSERT_EQUAL_TOL(gX,xxx,1e-8);
	
#endif
	
	return ctrans_mult(cx, cz);

}

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

