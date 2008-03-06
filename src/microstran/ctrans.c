
#include "ctrans.h"
#include <math.h>

#define PI 3.14159265358

ctrans_matrix ctrans_identity(void){
	ctrans_matrix c;
	unsigned i,j;
	for(i=0;i<4; ++i){
		for(j=0;j<4;++j){
			c.m[i][j] = 0.0;
		}
	};
	for(i=0;i<4;++i){
		c.m[i][i]=1.0;
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

ctrans_matrix ctrans_rotation(vec3 axis, double theta){
	double C = cos(theta);
	double S = sin(theta);
#define SQ(X) ((X)*(X))
	ctrans_matrix c;
	c = ctrans_identity();
	fprintf(stderr,"I =");
	ctrans_print(stderr,&c);
	// SOURCE: http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/index.htm
	c.m[0][0] += (1-C)*(SQ(axis.x)-1);
	c.m[0][1] = -axis.z*S+(1-C)*axis.x*axis.y;
	c.m[0][2] = axis.y*S+(1-C)*axis.x*axis.z;
	c.m[1][0] = axis.z*S+(1-C)*axis.x*axis.y;
	c.m[1][1] += (1-C)*(SQ(axis.y)-1);
	c.m[1][2] = -axis.x*S+(1-C)*axis.y*axis.z;
	c.m[2][0] = -axis.y*S+(1-C)*axis.x*axis.z;
	c.m[2][1] = axis.x*S+(1-C)*axis.y*axis.z;
	c.m[2][2] += (1-C)*(SQ(axis.z)-1);
	fprintf(stderr,"ROT =");
	ctrans_print(stderr,&c);
#undef SQ
	return c;
}

ctrans_matrix ctrans_rotation_axes(vec3 Z, vec3 X){
	ctrans_matrix cz, cx;
	vec3 z = vec3_create(0,0,1);
	vec3 dirz;
	double thetaz = vec3_angle_cross(z,Z,&dirz);
	fprintf(stderr,"Z-Rotation by %f around %f,%f,%f\n",thetaz*180./PI,dirz.x,dirz.y,dirz.z);
	cz = ctrans_rotation(dirz,thetaz);
	vec3 xd = vec3_rotate(X,dirz,thetaz);
	vec3 dirx;
	double thetax = vec3_angle_cross(xd,X,&dirx);
	if(thetax<1e-5){
		fprintf(stderr,"No X rotation required\n");
		return cz;
	}
	fprintf(stderr,"X-Rotation by %f around %f,%f,%f\n",thetax*180./PI,dirx.x,dirx.y,dirx.z);
	cx = ctrans_rotation(dirz,thetax);
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

