#include <string.h>
#include <assert.h>

#define MSTRANP_BUILD
#include "sections.h"
#include "new.h"
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
	PT(-bf/2, d/2);
	PT(-bf/2, d/2 - tf);
	PT(-tw/2, d/2 - tf);
	PT(-tw/2, -(d/2-tf));
	PT(-bf/2, -(d/2-tf));
	PT(-bf/2, -d/2);
	PT(bf/2, -d/2);
	PT(bf/2, -(d/2-tf));
	PT(tw/2, -(d/2-tf));
	PT(tw/2, d/2 - tf);
	PT(bf/2, d/2 - tf);
	PT(bf/2, d/2);
#undef PT

	for(i=0;i<12;++i)array_set(&(o->trace),i,&i);
	i = 0; // join back to the start again
	array_set(&(o->trace),12,&i);

	return o;
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
	}
	return n;
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


