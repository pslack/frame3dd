#include <microstran/ctrans.h>
#include "test.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358

class TestCtrans : public CppUnit::TestFixture{

public:

	void testidentity(){
		ctrans_matrix i = ctrans_identity();
		CPPUNIT_ASSERT(i.m[0][0]==1);
		CPPUNIT_ASSERT(i.m[1][0]==0);
		CPPUNIT_ASSERT(i.m[3][3]==1);
		CPPUNIT_ASSERT(i.m[2][3]==0);
		CPPUNIT_ASSERT(i.m[3][2]==0);
	}

	void testinv1(){
		ctrans_matrix i = ctrans_identity();
		ctrans_matrix ii = ctrans_inverse(i);
		CPPUNIT_ASSERT(ctrans_equal_tol(&i,&ii,1e-20));
	}

	void testinv2(){
		ctrans_matrix I = ctrans_identity();
		ctrans_matrix A = ctrans_identity();
		A.m[0][3] = 3.4; // add some off-diag values
		A.m[2][1] = 8.3;
		ctrans_matrix Ai = ctrans_inverse(A);
		ctrans_matrix AAi = ctrans_mult(A,Ai);
		CPPUNIT_ASSERT(ctrans_equal_tol(&AAi,&I,1e-20));

		ctrans_matrix AiA = ctrans_mult(Ai,A);
		CPPUNIT_ASSERT(ctrans_equal_tol(&AiA,&I,1e-20));
	}

	void testrotation(){
		vec3 axis = vec3_create(3,4,5);
		ctrans_matrix R = ctrans_rotation(axis,3.*PI/4);
		ctrans_matrix Ri = ctrans_rotation(axis,-3.*PI/4);
		ctrans_matrix RiR = ctrans_mult(Ri,R);
		ctrans_matrix I = ctrans_identity();
		CPPUNIT_ASSERT(ctrans_equal_tol(&RiR,&I,1e-15));
	}

	void testrot2(){
		vec3 axis = vec3_create(-3,-4,-5);
		ctrans_matrix R = ctrans_rotation(axis,7.*PI/4);
		ctrans_matrix Ri = ctrans_rotation(axis,-7.*PI/4);
		ctrans_matrix RiR = ctrans_mult(Ri,R);
		ctrans_matrix I = ctrans_identity();
		CPPUNIT_ASSERT(ctrans_equal_tol(&RiR,&I,1e-15));
	}

	void testrotz(){
		vec3 z = vec3_create(0,0,1);
		ctrans_matrix Rz = ctrans_rotation_z(-3.5*PI/4);
		ctrans_matrix R = ctrans_rotation(z, -3.5*PI/4);
		ctrans_matrix Ri = ctrans_inverse(R);
		ctrans_matrix RiRz = ctrans_mult(Ri,Rz);
		ctrans_matrix I = ctrans_identity();
		CPPUNIT_ASSERT(ctrans_equal_tol(&RiRz,&I,1e-15));
	}

	/// check consistency of vec3_rotate and ctrans_rotation for many random vectors
	void testrotrandom(){
		unsigned i,j;
		double maxlen = 30;
		for(i=0;i<100;++i){
			//fprintf(stderr,"i=%u...\n",i);
			vec3 axis = vec3_norm(vec3_create(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX));
			double theta = 2.*PI*rand()/RAND_MAX - PI;
			for(j=0;j<100;++j){
				// random vector to be rotated
				vec3 v = vec3_scale(
					vec3_norm(vec3_create(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX))
					,maxlen*float(rand())/RAND_MAX
				);
				CPPUNIT_ASSERT(vec3_equal_tol(ctrans_apply(ctrans_rotation(axis,theta),v),vec3_rotate(v,axis,theta),1e-13));
			}
		}
	}		
	
	/// test null transformation using ctrans_rotation_axes
	void testrotaxes(){
		vec3 Z = vec3_create(0,0,1);
		vec3 X = vec3_create(1,0,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		//CTRANS_PR(R);
		ctrans_matrix I = ctrans_identity();
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&I,1e-15));
	}

	/// test 180° rotation around Z axis using ctrans_rotation_axes
	void testrotaxes2(){
		vec3 Z = vec3_create(0,0,1);
		vec3 X = vec3_create(-1,0,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		//CTRANS_PR(R);
		
		/* equivalent rotation is a X and Y flip, Z unchanged */
		ctrans_matrix Re = ctrans_identity();
		Re.m[0][0] = -1;
		Re.m[1][1] = -1;

		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&Re,1e-15));
	}

	/// test 45° rotation around the +Z axis
	void testrotaxes3(){
		vec3 Z = vec3_create(0,0,1);
		vec3 X = vec3_norm(vec3_create(1,1,0));
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		ctrans_matrix R1 = ctrans_rotation(Z,PI/4);
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&R1,1e-8));
	}

	/// test 135° rotation around Z axis using ctrans_rotation_axes
	void testrotaxes4(){
		vec3 Z = vec3_create(0,0,1);
		vec3 X = vec3_norm(vec3_create(-1,1,0));
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		ctrans_matrix R1 = ctrans_rotation(Z,3*PI/4);
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&R1,1e-8));
	}

	/// test 45° rotation around X axis using ctrans_rotation_axes
	void testrotaxes5(){
		vec3 Z = vec3_norm(vec3_create(0,1,1));
		vec3 X = vec3_create(1,0,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		ctrans_matrix R1 = ctrans_rotation(X,-PI/4);
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&R1,1e-8));
	}

	/// test 135° rotation around X axis using ctrans_rotation_axes
	void testrotaxes6(){
		vec3 Z = vec3_norm(vec3_create(0,1,-1));
		vec3 X = vec3_create(1,0,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		ctrans_matrix R1 = ctrans_rotation(X,-3.*PI/4);
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&R1,1e-8));
	}

	void testrotaxes7(){
		vec3 X = vec3_create(-0.000000000000000000000000000000e+00, 9.908472411173366856118605028314e-01, 1.349879430547868131018418580425e-01);
		vec3 Z = vec3_create(1.000000000000000000000000000000e+00, 0.000000000000000000000000000000e+00, 0.000000000000000000000000000000e+00);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		CPPUNIT_ASSERT(vec3_equal_tol(ctrans_apply(R,vec3_create(1,0,0)),X,1e-8));
		CPPUNIT_ASSERT(vec3_equal_tol(ctrans_apply(R,vec3_create(0,0,1)),Z,1e-8));
	}		


#if 0
	// THIS TEST IS INVALID
	void testrotaxes7(){
		vec3 Z = vec3_create(1,0,0);
		vec3 X = vec3_create(0,1,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		CTRANS_PR(R);

		ctrans_matrix R1 = ctrans_rotation(vec3_create(1,1,1),-PI/3);
		CTRANS_PR(R1);

		fprintf(stderr,"Test ctrans_rotation_axes...\n");
		vec3 v1 = ctrans_apply(R,vec3_create(1,1,0)/* in X Y Z coors */);
		VEC3_PR(v1);
		CPPUNIT_ASSERT(vec3_equal_tol(v1,vec3_create(0,1,1)/* in x y z coors */,1e-8));


		fprintf(stderr,"Test ctrans_rotation about 1,1,1 by -PI/3...\n");
		vec3 unrot = ctrans_apply(R1,vec3_create(1,1,0));
		VEC3_PR(unrot);
		CPPUNIT_ASSERT(vec3_equal_tol(unrot,vec3_create(1,1,0),1e-8));

		CTRANS_PR(R1);
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&R1,1e-8));
	}		
#endif

	CPPUNIT_TEST_SUITE(TestCtrans);
	CPPUNIT_TEST(testidentity);
	CPPUNIT_TEST(testinv1);
	CPPUNIT_TEST(testinv2);
	CPPUNIT_TEST(testrotation);
	CPPUNIT_TEST(testrot2);
	CPPUNIT_TEST(testrotz);
	CPPUNIT_TEST(testrotrandom);
	CPPUNIT_TEST(testrotaxes);
	CPPUNIT_TEST(testrotaxes2);
	CPPUNIT_TEST(testrotaxes3);
	CPPUNIT_TEST(testrotaxes4);
	CPPUNIT_TEST(testrotaxes5);
	CPPUNIT_TEST(testrotaxes6);
	CPPUNIT_TEST_SUITE_END();

};

CPPUNIT_TEST_SUITE_REGISTRATION(TestCtrans);
