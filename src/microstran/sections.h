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
	,SECTION_ISEC
	,SECTION_SHS
	,SECTION_TOPHAT
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

/**
	I-sections (universal beams, universal columns, maybe other?)
*/
struct section_isec_struct{
	double m_linear; /* mass per linear metre */
	double d; /* section depth */
	double bf, tf; /* flange width and thickness */
	double tw; /* web thickness */
	double r1; /* root radius */
	double Ag; /* cross-section area */
	double lx, Zx, Sx, rx;
	double ly, Zy, Sy, ry;
	double J, lw; /* torsion constant and warping constant */
};

/**
	SHS square hollow sections
*/
struct section_shs_struct{
	double d; /* section width and depth (equal in all cases) */
	double t; /* section thickness */
	double m_linear;
	double Ag;
	double Ix, Zx, Zn, Sx, rx;
	double J, C;
	double k_f, lambda_e;
	char compactness;
	double Ze;
};

/**
	Top-hat sections
*/
struct section_tophat_struct{
	double a; /**< width at top of top hat (in m) */
	double b; /**< overall height of top hat (in m) */
	double c; /**< width of (a single) flange (in m) */
	double d; /**< length of lip (in m) */
	double theta; /**< angle of sidewalls to vertical (in radians) */
	double phi; /**< angle of lip to vertical (in radians) */
	double t; /**< sheetmetal thickness (in m) */
	/* don't know what structural properties to define here yet */
};

#define SECTION_NAME_MAX 40

struct section_struct{
	char name[SECTION_NAME_MAX];
	section_type type;
	union{
		struct section_chs_struct chs;
		struct section_isec_struct isec;
		struct section_shs_struct shs;
		struct section_tophat_struct tophat;
	};
};
typedef struct section_struct section;

/*
	Structure to return section outline to other codes
*/
typedef struct vec2_struct{
	double x,y;
} vec2;

typedef struct section_outline_struct {
	array point; /**< (of vec2) vertices in the cross-section */
	array trace; /**< (of int) array of point IDs, or -1 to pick up 'pen' */
} section_outline;

typedef struct section_library_struct{
	array a;
} section_library;

MSTRANP_API void section_outline_destroy(section_outline *o);

MSTRANP_API section_library *section_library_create();
MSTRANP_API void section_library_destroy(section_library *l);

MSTRANP_API const section *section_find(const section_library *l, const char *name);

/* chs routines */
MSTRANP_API int section_is_chs(const section *s);
MSTRANP_API double section_chs_outside_diameter(const section *s);
MSTRANP_API double section_chs_thickness(const section *s);

/* i-section routines */
MSTRANP_API int section_is_isec(const section *s);
MSTRANP_API double section_isec_depth(const section *s);
MSTRANP_API double section_isec_flange_width(const section *s);
MSTRANP_API double section_isec_flange_thickness(const section *s);
MSTRANP_API double section_isec_web_thickness(const section *s);
MSTRANP_API section_outline *section_isec_outline(const section *s);

/* shs routines */
MSTRANP_API int section_is_shs(const section *s);
MSTRANP_API double section_shs_depth(const section *s);
MSTRANP_API double section_shs_thickness(const section *s);
MSTRANP_API section_outline *section_shs_outline(const section *s);

/* tophat routines */
MSTRANP_API int section_is_tophat(const section *s);
MSTRANP_API double section_tophat_top_width(const section *s);
MSTRANP_API double section_tophat_depth(const section *s);
MSTRANP_API double section_tophat_flange_width(const section *s);
MSTRANP_API double section_tophat_lip_width(const section *s);
MSTRANP_API double section_tophat_theta(const section *s);
MSTRANP_API double section_tophat_phi(const section *s);
MSTRANP_API double section_tophat_thickness(const section *s);
MSTRANP_API section_outline *section_tophat_outline(const section *s);

MSTRANP_API int section_print(FILE *f, const section *s);

#ifdef __cplusplus
}
#endif

#endif /* SECTIONS_H */
