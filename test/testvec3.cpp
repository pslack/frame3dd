#include <microstran/vec3.h>
#include "test.h"

#include <math.h>

#define PI 3.14159265358

class TestVec3 : public CppUnit::TestFixture{

public:

	void test1(){
		vec3 v;
		v = vec3_create(1,2,3);
		CPPUNIT_ASSERT(v.x==1);
		CPPUNIT_ASSERT(v.y==2);
		CPPUNIT_ASSERT(v.z==3);
	}

	void testadd(){
		vec3 v1 = vec3_create(1.,2.,3.);
		vec3 v2 = vec3_create(4.,5.,6.);
		CPPUNIT_ASSERT(vec3_equal(vec3_add(v1,v2), vec3_create(5,7,9)));
	}

	void testdiff(){
		vec3 v1 = vec3_create(5.,-1.,6.);
		vec3 v2 = vec3_create(2.,4.,-8.);
		CPPUNIT_ASSERT(vec3_equal(vec3_diff(v1,v2),vec3_create(3.,-5.,14.)));
	}

	void testdistance1(){
		vec3 v1 = vec3_create(1,1,0);
		vec3 v2 = vec3_create(0,0,0);
		CPPUNIT_ASSERT_EQUAL(vec3_mod(vec3_diff(v1,v2)),sqrt(2));
	}

	void testdistance2(){
		vec3 v1 = vec3_create(1,-1,0);
		vec3 v2 = vec3_create(0,0,0);
		CPPUNIT_ASSERT_EQUAL(vec3_mod(vec3_diff(v1,v2)),sqrt(2));
	}

	void testdistance3(){
		vec3 v1 = vec3_create(3,3,3);
		vec3 v2 = vec3_create(0,0,0);
		CPPUNIT_ASSERT_EQUAL(vec3_mod(vec3_diff(v1,v2)),sqrt(27));
	}

	void testdot(){
		vec3 v1 = vec3_create(1,2,3);
		vec3 v2 = vec3_create(4,5,6);
		double d = vec3_dot(v1,v2);
		CPPUNIT_ASSERT(vec3_equal(v1,vec3_create(1,2,3)));
		CPPUNIT_ASSERT(vec3_equal(v2,vec3_create(4,5,6)));
		CPPUNIT_ASSERT(d==4+10+18);
	}

	void testcross(){
		vec3 x = vec3_create(1,0,0);
		vec3 y = vec3_create(0,1,0);
		CPPUNIT_ASSERT(vec3_equal(vec3_cross(x,y),vec3_create(0,0,1)));

		vec3 z = vec3_create(0,0,1);
		CPPUNIT_ASSERT(vec3_equal(vec3_cross(x,z),vec3_scale(y,-1)));

		/// @TODO more tests here
	}

	void testangle(){
		vec3 x = vec3_create(1,0,0);
		vec3 y = vec3_create(0,1,0);
		CPPUNIT_ASSERT(fabs( vec3_angle(x,y) - PI/2. ) < 1e-8 );

		CPPUNIT_ASSERT(fabs( vec3_angle(vec3_add(x,y),y) - PI/4. ) < 1e-8);
	}

	/**
		A few tests of the behaviour of cross product/angle functions with
		large/reflex angles
	*/
	void testcrossangle(){
		//  y
		//   \._ x
		//
		vec3 x = vec3_create(1,0,0);
		vec3 y = vec3_create(-1,1,0);
		CPPUNIT_ASSERT(fabs(vec3_angle(x,y) - 3./4*PI) < 1e-8 );
		vec3 c = vec3_cross(x,y);
		CPPUNIT_ASSERT(fabs( c.z - 1 ) < 1e-8 );

		//
		//     ._ x
		// y2 /

		vec3 y2 = vec3_create(-1,-1,0);
		CPPUNIT_ASSERT(fabs(vec3_angle(x,y2) - 3./4*PI) < 1e-8 );
		vec3 c2 = vec3_cross(x,y2);
		CPPUNIT_ASSERT(fabs( c2.z + 1 ) < 1e-8 );

		//
		//     ._ x
		// y2 /

		vec3 c3 = vec3_cross(y2,x);
		CPPUNIT_ASSERT(fabs(vec3_angle(y2,x) - 3./4*PI) < 1e-8 ); // commutative
		CPPUNIT_ASSERT(fabs( c3.z - 1 ) < 1e-8 );

	}

#if 0
	/// @TODO need more tests for the 'rotate' routine!
	void testrotate(){
		vec3 x = vec3_create(1,0,0);
		vec3 y = vec3_create(0,1,0);
		VEC3_PR(y);
		y = vec3_rotate(y,x,PI/2);
		VEC3_PR(y);
		CPPUNIT_ASSERT(vec3_equal_tol(y,vec3_create(0.,0.,1.),1e-11));

		vec3 r = vec3_create(1,1,1);
		vec3 z = vec3_create(0,0,1);
		z = vec3_rotate(z,r,120.*PI/180);
		CPPUNIT_ASSERT(vec3_equal(z,vec3_create(1,0,0)));
	}
#endif

#if 0
	void testpolar(){
		Vector x(1,0,0);
		Vector r = Vector::polar(1,45*DEG,90*DEG);
		x.rotate(Vector(0,0,1),45);
		CPPUNIT_ASSERT(r==x);
	}
#endif

	CPPUNIT_TEST_SUITE(TestVec3);
	CPPUNIT_TEST(test1);
	CPPUNIT_TEST(testadd);
	CPPUNIT_TEST(testdiff);
	CPPUNIT_TEST(testdistance1);
	CPPUNIT_TEST(testdistance2);
	CPPUNIT_TEST(testdistance3);
	CPPUNIT_TEST(testdot);
	CPPUNIT_TEST(testcross);
	CPPUNIT_TEST(testangle);
	CPPUNIT_TEST(testcrossangle);
	CPPUNIT_TEST_SUITE_END();

};

CPPUNIT_TEST_SUITE_REGISTRATION(TestVec3);

