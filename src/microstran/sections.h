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
	,SECTION_ROD
} section_type;

/**
	CHS section properties
*/
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
	char inverted; /**< direction of section in ±Y direction 0=open downwards, 1= open upwards */
	/* don't know what structural properties to define here yet */
};

/**
	CHS section properties
*/
struct section_rod_struct{
	double d; /**< diameter (in m) */
	double A; /**< area (m²) */
	double Ax;
	double J;
	double Ix;
	double M;
};

#define SECTION_NAME_MAX 40

/**
	General section properties structure
*/
struct section_struct{
	char name[SECTION_NAME_MAX];
	section_type type;
	union{
		struct section_chs_struct chs;
		struct section_isec_struct isec;
		struct section_shs_struct shs;
		struct section_tophat_struct tophat;
		struct section_rod_struct rod;
	};
};
typedef struct section_struct section;

/**
	Structure to return section outline to other codes
*/
typedef struct vec2_struct{
	double x,y;
} vec2;

/**
	Section outline (trace of points for use in rendering steel sections)

	The section is defined as a number of vertices (of type vec2) which are
	indexed in the array 'point'. The section outline is then created as a
	series of these indices. For each edge that is traced, the 'stuff' is 
	assumed to be on the RIGHT of the line in the direction of travel.

	If there are discontinuities in the section, these are represented by a
	'pen up' code (=-1) in the trace array.
*/
typedef struct section_outline_struct {
	array point; /**< (of vec2) vertices in the cross-section */
	array trace; /**< (of int) array of point IDs, or -1 to pick up 'pen' */
} section_outline;

/**
	Array of known steel sections of different sizes and types, as loaded
	using section_library_create and related parser functions.
*/
typedef struct section_library_struct{
	array a;
} section_library;

MSTRANP_API void section_outline_destroy(section_outline *o);

/**
	Create a copy of a section_outline with inverted y coordinates.
	This function also re-orders the 'trace' array, to obey the 'stuff-to-the-
	right' rule for the section outline, required for accurate graphical 
	rendering.
*/
MSTRANP_API section_outline  *section_outline_copy_inverted(section_outline *o);

/**
	Write the data contained in a section outline for debuggin purposes.
*/
MSTRANP_API int section_outline_print(FILE *f, const section_outline *o);

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
MSTRANP_API char section_tophat_inverted(const section *s);
MSTRANP_API section_outline *section_tophat_outline(const section *s);

/* chs routines */
MSTRANP_API int section_is_rod(const section *s);
MSTRANP_API double section_rod_diameter(const section *s);

MSTRANP_API int section_print(FILE *f, const section *s);

/**
	Get section outline, if available for this section.
	If not available (eg for CHS), return NULL (in which case other
	methods need to be used for graphical rendering)
*/
MSTRANP_API section_outline *section_get_outline(const section *s);

/** approximate diameter of section for use in fast rendering routines */
MSTRANP_API double section_outline_approx_diameter(const section_outline *o);

#ifdef __cplusplus
}
#endif

#endif /* SECTIONS_H */
