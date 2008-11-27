/*
 FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://www.duke.edu/~hpgavin/frame/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2008  Henri P. Gavin
 
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

#ifndef FRAME_COORDTRANS_H
#define FRAME_COORDTRANS_H

#include "common.h"

/**
	COORD_TRANS -  evaluate the 3D coordinate transformation coefficients 1dec04
	Default order of coordinate rotations...  typical for Y as the vertical axis
	1. rotate about the global Z axis
	2. rotate about the global Y axis
	3. rotate about the local  x axis --- element 'roll'

	If Zvert is defined as 1, then the order of coordinate rotations is typical
	for Z as the vertical axis
	1. rotate about the global Y axis
	2. rotate about the global Z axis
	3. rotate about the local  x axis --- element 'roll'

	Q=TF;   U=TD;   T'T=I;   Q=kU;   TF=kTD;   T'TF=T'kTD;   T'kT = K;   F=KD
*/
void coord_trans(
	float *x, float *y, float *z
	, float L
	, int j1, int j2
	, float *t1, float *t2, float *t3, float *t4, float *t5
	, float *t6, float *t7, float *t8, float *t9
	, float p /**< the roll angle (radians) */
);

#endif /* FRAME_COORDTRANS_H */

