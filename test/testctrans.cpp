#include <microstran/ctrans.h>
#include "test.h"

#include <math.h>

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
	
	void testrotaxes(){
		vec3 Z = vec3_create(0,0,1);
		vec3 X = vec3_create(1,0,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		CTRANS_PR(R);
		ctrans_matrix I = ctrans_identity();
		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&I,1e-15));
	}

	void testrotaxes2(){
		vec3 Z = vec3_create(0,0,1);
		vec3 X = vec3_create(-1,0,0);
		ctrans_matrix R = ctrans_rotation_axes(Z,X);
		CTRANS_PR(R);
		
		/* equivalent rotation is a X and Y flip, Z unchanged */
		ctrans_matrix Re = ctrans_identity();
		Re.m[0][0] = -1;
		Re.m[1][1] = -1;

		CPPUNIT_ASSERT(ctrans_equal_tol(&R,&Re,1e-15));
	}

	CPPUNIT_TEST_SUITE(TestCtrans);
	CPPUNIT_TEST(testidentity);
	CPPUNIT_TEST(testinv1);
	CPPUNIT_TEST(testinv2);
	CPPUNIT_TEST(testrotation);
	CPPUNIT_TEST(testrot2);
	CPPUNIT_TEST(testrotz);
	CPPUNIT_TEST(testrotaxes);
	CPPUNIT_TEST(testrotaxes2);
	CPPUNIT_TEST_SUITE_END();

};

CPPUNIT_TEST_SUITE_REGISTRATION(TestCtrans);
