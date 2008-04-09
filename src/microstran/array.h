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
*//** @file
	data type for dynamically-size array in C. this array will
	automatically resize when you use array_set outside its bounds. it is
	dense, however: storage is allocated for members in the range 0..num.

	you can't shrink the array unless you just copy a range into a new array
	then destroy the old one.

	note that when the array is resized, data is copied to a new location, so
	pointers into the array are a bad idea. you should used array indices
	instead.
*/

#ifndef MSTRANP_ARRAY_H
#define MSTRANP_ARRAY_H

#include "config.h"
#include <malloc.h>

#ifdef __cplusplus
extern "C"{
#endif

typedef struct array_{
	size_t eachsize;
	unsigned num;
	unsigned cap;
	void *data;
} array;

/** create a new array. prefer to use ARRAY_CREATE instead. */
array array_create(size_t eachsize, unsigned cap);

#define ARRAY_CREATE(TYPE,CAP) array_create(sizeof(TYPE),(CAP))

#define ARRAY_NUM(ARRAY) ((ARRAY).num)

/** set a value in the array. pass value by pointer; it is *copied* into place. @return pointer to the updated data element */
void *array_set(array *a, unsigned index, void *val);

/** look up a value in the array. @return pointer to the location where the data is stored. */
MSTRANP_API void *array_get(array *a, unsigned index);

/** place a value at the end of the array. value is copied into place. @return a pointer to the new member */
void *array_append(array *a, void *val);

/** destroy memory contained in array and clear num/cap of the array struct */
void array_destroy(array *a);

/** copy an array. capacity of new array will include no spare space */
array array_copy(const array a);

#ifdef __cplusplus
};
#endif

#endif /* MSTRANP_ARRAY_H */

