#include "coordtrans.h"

void coord_trans(
		float *x, float *y, float *z
		, float L
		, int j1, int j2
		, float *t1, float *t2, float *t3, float *t4, float *t5
		, float *t6, float *t7, float *t8, float *t9
		, float p /**< the roll angle (radians) */
){
	float	Cx, Cy, Cz, den,		/* direction cosines	*/
		Cp, Sp;			/* cosine and sine of roll angle */

	Cx = (x[j2] - x[j1]) / L;
	Cy = (y[j2] - y[j1]) / L;
	Cz = (z[j2] - z[j1]) / L;

	*t1 = *t2 = *t3 = *t4 = *t5 = *t6 = *t7 = *t8 = *t9 = 0.0;

	Cp = cos(p);
	Sp = sin(p);

#if Zvert

	if ( fabs(Cz) == 1.0 ) {
		*t3 =  Cz;
		*t4 = -Cz*Sp;
		*t5 =  Cp;
		*t7 = -Cz*Cp;
		*t8 = -Sp;
	} else {

		den = sqrt ( 1.0 - Cz*Cz );

		*t1 = Cx;
	   	*t2 = Cy;
		*t3 = Cz;

		*t4 = (-Cx*Cz*Sp - Cy*Cp)/den;    
 		*t5 = (-Cy*Cz*Sp + Cx*Cp)/den;
		*t6 = Sp*den;

		*t7 = (-Cx*Cz*Cp + Cy*Sp)/den;
		*t8 = (-Cy*Cz*Cp - Cx*Sp)/den;
	   	*t9 = Cp*den;
	}

#else

	if ( fabs(Cy) == 1.0 ) {
		*t2 =  Cy;
		*t4 = -Cy*Cp;
		*t6 =  Sp;
		*t7 =  Cy*Sp;
		*t9 =  Cp;
	} else {

		den = sqrt ( 1.0 - Cy*Cy );

		*t1 = Cx;
	   	*t2 = Cy;
		*t3 = Cz;

		*t4 = (-Cx*Cy*Cp - Cz*Sp)/den;    
		*t5 = den*Cp;
 		*t6 = (-Cy*Cz*Cp + Cx*Sp)/den;

		*t7 = (Cx*Cy*Sp - Cz*Cp)/den;
	   	*t8 = -den*Sp;
		*t9 = (Cy*Cz*Sp + Cx*Cp)/den;
	}

#endif

	return;
}

