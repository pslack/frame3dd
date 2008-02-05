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


