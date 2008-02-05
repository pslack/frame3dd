/** @FILE
	routines to load a section library from a text file such as
	'properties.txt' provided here.
*/

/* define all the 'known' section types */

#ifndef SECTIONS_H
#define SECTIONS_H

#ifdef __cplusplus
extern "C"{
#endif

#include "config.h"
#include "array.h"
#include <stdio.h>

typedef enum{
	SECTION_CHS=0
} section_type;

struct section_chs_struct{
	double d_o;
	double t;
	double A;
	double I;
	double Z;
	double S;
	double r;
	double J;
	double C;
	double k_f;
	double lambda_s;
	char compactness;
	double Ze;
};

#define SECTION_NAME_MAX 40

struct section_struct{
	char name[SECTION_NAME_MAX];
	section_type type;
	union{
		struct section_chs_struct chs;
	};
};

typedef struct section_struct section;

typedef struct section_library_struct{
	array a;
} section_library;

MSTRANP_API section_library *section_library_create();
MSTRANP_API void section_library_destroy(section_library *l);

MSTRANP_API const section *section_find(const section_library *l, const char *name);

MSTRANP_API int section_is_chs(const section *s);
MSTRANP_API double section_chs_outside_diameter(const section *s);
MSTRANP_API double section_chs_thickness(const section *s);

MSTRANP_API int section_print(FILE *f, const section *s);

#ifdef __cplusplus
}
#endif

#endif /* SECTIONS_H */
