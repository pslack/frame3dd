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
*/
#include <string.h>
#include <assert.h>
#include <math.h>

#define MSTRANP_BUILD
#include "sections.h"
#include "new.h"
#include "array.h"

#define PI	3.141592653589793

/**
	Look up a section from the properties library by name
	@param library an array of available sections to search throu
	@param name name of the section to look for
	@return pointer to the section found, or NULL if not found.
*/
const section *section_find(const section_library *library, const char *name){
	unsigned i;
	unsigned n = ARRAY_NUM(library->a);
	const section *s;
	for(i=0; i<n; ++i){
		s = (const section *)array_get((array *)&(library->a),i);
		if(strcmp(s->name,name)==0){
			return s;
		}
	}
	return NULL;
}

/* CHS routines */

int section_is_chs(const section *s){
	if(s->type==SECTION_CHS)return 1;
	return 0;
}

double section_chs_outside_diameter(const section *s){
	assert(section_is_chs(s));
	return s->chs.d_o;
}

double section_chs_thickness(const section *s){
	assert(section_is_chs(s));
	return s->chs.t;
}

/* I-SECTION routines */

int section_is_isec(const section *s){
	if(s->type==SECTION_ISEC)return 1;
	return 0;
}

double section_isec_depth(const section *s){
	assert(section_is_isec(s));
	return s->isec.d;
}

double section_isec_flange_width(const section *s){
	assert(section_is_isec(s));
	return s->isec.bf;
}

double section_isec_flange_thickness(const section *s){
	assert(section_is_isec(s));
	return s->isec.tf;
}

double section_isec_web_thickness(const section *s){
	assert(section_is_isec(s));
	return s->isec.tw;
}

/* SHS routines */

int section_is_shs(const section *s){
	if(s->type==SECTION_SHS)return 1;
	return 0;
}

double section_shs_depth(const section *s){
	assert(section_is_shs(s));
	return s->shs.d;
}

double section_shs_thickness(const section *s){
	assert(section_is_shs(s));
	return s->shs.t;
}

/* Tophat routines */

int section_is_tophat(const section *s){
	if(s->type==SECTION_TOPHAT)return 1;
	return 0;
}
	
double section_tophat_top_width(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.a;
}

double section_tophat_depth(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.b;
}

double section_tophat_flange_width(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.c;
}

double section_tophat_lip_width(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.d;
}

double section_tophat_theta(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.theta;
}

double section_tophat_phi(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.phi;
}

char section_tophat_inverted(const section *s){
	assert(section_is_tophat(s));
	return s->tophat.inverted;
}

/* Rod routines */

int section_is_rod(const section *s){
	if(s->type==SECTION_ROD)return 1;
	return 0;
}

double section_rod_diameter(const section *s){
	assert(section_is_rod(s));
	return s->rod.d;
}

/* utility routine for use with array_set */
static vec2 *vec2_set(vec2 *v, double x, double y){
	v->x = x;
	v->y = y;
	return v;
}

/**
	Output a data structure containing the perimeter of the steel section.
	We can use this to generate 3D 'prism' geometry elsewhere.
*/
section_outline *section_isec_outline(const section *s){
	section_outline *o;
	vec2 v;
	assert(section_is_isec(s));
	double tf = s->isec.tf;
	double tw = s->isec.tw;
	double bf = s->isec.bf;
	double d = s->isec.d;
	
	o = NEW(section_outline);	
	o->point = ARRAY_CREATE(vec2,12);
	o->trace = ARRAY_CREATE(int,12);
	int i=0;
#define PT(X,Y) array_set(&(o->point),i++,vec2_set(&v,(X),(Y)))
	/* trace around the contour in a clockwise direction, following the
	convention that the surface 'exterior' is on the left side of the contour */
	PT(bf/2, d/2); /* 0 */
	PT(bf/2, d/2 - tf); /* 1 */
	PT(tw/2, d/2 - tf); /* 2 */
	PT(tw/2, -(d/2-tf)); /* 3 */
	PT(bf/2, -(d/2-tf)); /* 4 */
	PT(bf/2, -d/2); /* 5 */
	PT(-bf/2, -d/2); /* 6 */
	PT(-bf/2, -(d/2-tf)); /* 7 */
	PT(-tw/2, -(d/2-tf)); /* 8 */
	PT(-tw/2, d/2 - tf); /* 9 */
	PT(-bf/2, d/2 - tf); /* 10 */
	PT(-bf/2, d/2); /* 11 */
#undef PT

	for(i=0;i<12;++i)array_set(&(o->trace),i,&i);
	i = 0; // join back to the start again
	array_set(&(o->trace),12,&i);

	return o;
}


/**
	Output a data structure containing the perimeter of the steel section.
	We can use this to generate 3D 'prism' geometry elsewhere.
*/
section_outline *section_shs_outline(const section *s){
	section_outline *o;
	vec2 v;
	assert(section_is_shs(s));
	double d = s->shs.d;
	
	o = NEW(section_outline);	
	o->point = ARRAY_CREATE(vec2,4);
	o->trace = ARRAY_CREATE(int,4);
	int i=0;
#define PT(X,Y) array_set(&(o->point),i++,vec2_set(&v,(X),(Y)))
	/* trace around the contour in a clockwise direction, following the
	convention that the surface 'exterior' is on the left side of the contour */
	PT(d/2, d/2); /* 0 */
	PT(d/2, -d/2); /* 1 */
	PT(-d/2, -d/2); /* 2 */
	PT(-d/2, d/2); /* 3 */
#undef PT

	for(i=0;i<4;++i)array_set(&(o->trace),i,&i);
	i = 0; // join back to the start again
	array_set(&(o->trace),4,&i);
	return o;
}

/**
	Output a data structure containing the perimeter of the TOP HAT section.
	We can use this to generate 3D 'prism' geometry elsewhere.
*/
section_outline *section_tophat_outline(const section *s){
	section_outline *o, *oi;
	vec2 v;
	assert(section_is_tophat(s));
	double a = s->tophat.a;
	double b = s->tophat.b;
	double c = s->tophat.c;
	double d = s->tophat.d;
	double t = s->tophat.t;
	double theta = s->tophat.theta;
	double phi = s->tophat.phi;

	double bx = (b - t)*tan(theta);
	double c1 = c - t/tan(0.5*theta + 0.25*PI) - t/tan(0.5*phi + 0.25*PI);

	double dx = d*cos(phi);
	double dx1 = dx - t / tan(0.5*phi + 0.25*PI);
	
	o = NEW(section_outline);	
	o->point = ARRAY_CREATE(vec2,16);
	o->trace = ARRAY_CREATE(int,17);
	int i=0;
#define PT(X,Y) array_set(&(o->point),i++,vec2_set(&v,(X),(Y)))
	/* trace around the contour in a clockwise direction, following the
	convention that the surface 'exterior' is on the left side of the contour */
	/* these points are plotted relative to an origin at the middle of the plane
	containing the flanges */
	PT(a/2, b); /* 0 */
	PT(a/2 + bx, t); /* 1 */
	PT(a/2 + bx + c1, t); /* 2 */
	PT(a/2 + bx + c1 + dx1*sin(phi), t + dx1*cos(phi)); /* 3 */
	PT(a/2 + bx + c1 + dx1*sin(phi) + t * cos(phi), t + dx1*cos(phi) - t*sin(phi)); /* 4 */
	PT(a/2 + bx - t/tan(0.5*theta + 0.25*PI) + c, 0); /* 5 */
	PT(a/2 + bx - t/tan(0.5*theta + 0.25*PI), 0); /* 6 */
	PT(a/2 - t/tan(0.5*theta + 0.25*PI), b - t); /* 7 */

	PT(-(a/2 - t/tan(0.5*theta + 0.25*PI)), b - t); /* 7 */
	PT(-(a/2 + bx - t/tan(0.5*theta + 0.25*PI)), 0); /* 6 */
	PT(-(a/2 + bx - t/tan(0.5*theta + 0.25*PI) + c), 0); /* 5 */
	PT(-(a/2 + bx + c1 + dx1*sin(phi) + t * cos(phi)), t + dx1*cos(phi) - t*sin(phi)); /* 4 */
	PT(-(a/2 + bx + c1 + dx1*sin(phi)), t + dx1*cos(phi)); /* 3 */
	PT(-(a/2 + bx + c1), t); /* 2 */
	PT(-(a/2 + bx), t); /* 1 */
	PT(-(a/2), b); /* 0 */
#undef PT

	for(i=0;i<16;++i)array_set(&(o->trace),i,&i);
	i = 0; // join back to the start again
	array_set(&(o->trace),16,&i);

	double cx, cy;
	cx = 0;
	cy = 0.5 * b;
	
	/* adjust for location of the centroid */
	unsigned n = ARRAY_NUM(o->point);
	for(i=0;i<n;++i){
		vec2 *p = (vec2 *)array_get(&(o->point),i);
		p->x -= cx;
		p->y -= cy;
	}

	/* invert if necessary */
	if(s->tophat.inverted){
		oi = section_outline_copy_inverted(o);
		section_outline_destroy(o);
		return oi;
	}

	return o;
}

section_outline *section_outline_copy_inverted(section_outline *o){
	/* to copy an inverted section outline, we must reverse the order of the
	 	points, because otherwise the section will also be inside out. */

	section_outline *o1;
	int *t;
	unsigned i, n = ARRAY_NUM(o->trace);
	vec2 *v;

	o1 = NEW(section_outline);
	o1->point = array_copy(o->point);
	o1->trace = ARRAY_CREATE(int,n);

	for(i=0; i<n; ++i){
		t = (int *)array_get(&o->trace, n - 1 - i);
		array_set(&(o1->trace),i,t);
	}

	// now flip the up-down coordinates
	n = ARRAY_NUM(o1->point);
	for(i=0; i<n; ++i){
		v = array_get(&o1->point, i);
		v->y = -v->y;
	}

	return o1;
}

void section_outline_destroy(section_outline *o){
	array_destroy(&(o->point));
	array_destroy(&(o->trace));
	free(o);
}

/* general routines */

int section_print(FILE *f, const section *s){
	int n = 0;
	if(section_is_chs(s)){
		n += fprintf(f,"%s\n",s->name);
		n += fprintf(f,"\tCHS, outside diameter = %f, thickness = %f\n", s->chs.d_o, s->chs.t);
	}else if(section_is_shs(s)){
		n += fprintf(f,"%s\n",s->name);
		n += fprintf(f,"\tSHS, width and depth = %f, thickness = %f\n", s->shs.d, s->shs.t);
	}else if(section_is_tophat(s)){
		n += fprintf(f,"%s\n",s->name);
		n += fprintf(f,"\tTOPHAT, a = %f, b = %f, c = %f, d = %f, thickness = %f%s\n"
			, s->tophat.a, s->tophat.b, s->tophat.c, s->tophat.d
			, s->tophat.t, (s->tophat.inverted? ", INVERTED" : "")
		);
	}else if(section_is_rod(s)){
		n += fprintf(f,"%s\n",s->name);
		n += fprintf(f,"\tROD, diameter = %f\n",s->rod.d);
	}
	return n;
}

#define SQ(X) ((X)*(X))

#if 0
double section_approx_radius(const section *s){
	switch(s->type){
		case SECTION_CHS:
			return s->chs.d_o;
		case SECTION_ISEC:
			return sqrt(SQ(s->isec.d)+SQ(s->isec.bf));
		case SECTION_SHS:
			return sqrt(2)*s->shs.d;
		case SECTION_TOPHAT:
			return sqrt(SQ(s->tophat.b)+SQ(s->tophat.a+s->tophat.c));
	}
}
#endif

double section_outline_approx_diameter(const section_outline *o){
	unsigned i, n = ARRAY_NUM(o->point);
	double maxx=0, maxy=0;
	for(i=0;i<n;++i){
		vec2 *p = (vec2 *)array_get((array *)&(o->point),i);
		if(fabs(p->x) > maxx)maxx=fabs(p->x);
		if(fabs(p->y) > maxy)maxy=fabs(p->y);
	}
	return 2.*sqrt(SQ(maxx)+SQ(maxy));
}

section_outline *section_get_outline(const section *s){
	switch(s->type){
		case SECTION_ISEC: return section_isec_outline(s);
		case SECTION_SHS: return section_shs_outline(s);
		case SECTION_TOPHAT: return section_tophat_outline(s);
		default:
			return NULL;
	}
}

int section_outline_print(FILE *f, const section_outline *o){
	unsigned i;
	int t;
	vec2 *v;
	unsigned n;
	unsigned c = 0;
	n = ARRAY_NUM(o->trace);
	c += fprintf(f,"%d trace points on %d vertices:\n",ARRAY_NUM(o->trace),ARRAY_NUM(o->point));
	int from=0;
	for(i=0; i<n; ++i){
		t = *(int *)array_get((array *)&o->trace,i);
		if(t==-1){
			from = 0; 
			c += fprintf(f," . . .\n");
			continue;
		}
		v = (vec2 *)array_get((array *)&o->point,t);
		if(from){
			c += fprintf(f," --> ");
		}else{
			c += fprintf(f,"     ");
		}

		c += fprintf(f,"%d:(%f,%f)\n", t, v->x, v->y);
		from = 1;
	}
	c += fprintf(f,"\n");
	return c;
}

section_library *section_library_create(){
	section_library *l;
	l = NEW(section_library);
	l->a = ARRAY_CREATE(section,20);
	return l;
}

void section_library_destroy(section_library *l){
	array_destroy(&(l->a));
	free(l);
}


