/*	
 FRAME3DD: 
 Static and dynamic structural analysis of 2D & 3D frames and trusses
 with elastic and geometric stiffness.
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
*//** @file
	this file contains some #defines to set up 'float' to be 'double' instead.
*/

#ifndef FRAME_COMMON_H
#define FRAME_COMMON_H


#ifndef VERSION
# define VERSION "20080909"
#endif


#ifndef PI
#define PI 3.14159265358979323846
#endif


#define Zvert 1		/* Zvert=0: Y axis is vertical ... rot Z, then rot Y */
			/* Zvert=1: Z axis is vertical ... rot Y, then rot Z */


#endif /* FRAME_COMMON_H */

