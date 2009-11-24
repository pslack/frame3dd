/*
 FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2009  Henri P. Gavin
 
    FRAME3DD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FRAME3DD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "coordtrans.h"

/* -------------------------------------------------------------------------
COORD_TRANS - calculate the 9 elements of the block-diagonal 12-by-12
coordinate transformation matrix, t1, t2, ..., t9.  

These coordinate transformation factors are used to:
* transform frame element end forces from the element (local) coordinate system 
to the structral (global) coordinate system.  
* transfrom end displacements from the structural (global) coordinate system
to the element (local) coordinate system,
* transform the frame element stiffness and mass matrices
from element (local) coordinates to structral (global) coordinates.

Element matrix coordinate transformations are carried out by function ATMA
in frame3dd.c

Currently coordinate transformations do not consider the effect of 
finite joint sizes ... this needs work, and could require a substantial
re-write of much of the code.  

Currently the effect of finite joint sizes is used only in the calculation
of the element stiffness matrices.  
------------------------------------------------------------------------- */
void coord_trans(
		vec3 *xyz,
		double L,
		int j1, int j2,
		double *t1, double *t2, double *t3, double *t4, double *t5,
		double *t6, double *t7, double *t8, double *t9,
		float p			/**< the roll angle (radians) */
){
	double	Cx, Cy, Cz, den,		/* direction cosines	*/
		Cp, Sp;			/* cosine and sine of roll angle */

	Cx = (xyz[j2].x - xyz[j1].x) / L;
	Cy = (xyz[j2].y - xyz[j1].y) / L;
	Cz = (xyz[j2].z - xyz[j1].z) / L;

	*t1 = *t2 = *t3 = *t4 = *t5 = *t6 = *t7 = *t8 = *t9 = 0.0;

	Cp = cos(p);
	Sp = sin(p);

#if Zvert				// the global Z axis is vertical

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

#else					// the global Y axis is vertical

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

