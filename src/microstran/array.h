/** @FILE
	data type for dynamically-size array in C. this array will
	automatically resize when you use array_set outside its bounds. it is
	dense, however: storage is allocated for members in the range 0..num.

	you can't shrink the array unless you just copy a range into a new array
	then destroy the old one.

	note that when the array is resized, data is copied to a new location, so
	pointers into the array are a bad idea. you should used array indices 
	instead.
*/
#ifndef ARRAY_H
#define ARRAY_H

#include <malloc.h>

typedef struct array_{
	size_t eachsize;
	unsigned num;
	unsigned cap;
	void *data;
} array;

/** create a new array. prefer to use ARRAY_CREATE instead. */
array array_create(size_t eachsize, unsigned cap);

#define ARRAY_CREATE(TYPE,CAP) array_create(sizeof(TYPE),(CAP))

#define ARRAY_NUM(ARRAY) (ARRAY.num)

/** set a value in the array. pass value by pointer; it is *copied* into place. @return pointer to the updated data element */
void *array_set(array *a, unsigned index, void *val);

/** look up a value in the array. @return pointer to the location where the data is stored. */
void *array_get(array *a, unsigned index);

/** place a value at the end of the array. value is copied into place. @return a pointer to the new member */
void *array_append(array *a, void *val);

/** destroy memory contained in array and clear num/cap of the array struct */
void array_destroy(array *a);

/** copy an array. capacity of new array will include no spare space */
array array_copy(const array a);

#endif

