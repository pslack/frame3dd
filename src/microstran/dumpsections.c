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
*//**@FILE
	little program to demonstrate that we can meaninfully parse the 'properties.txt' file
*/

#include "sections.h"
#include "sectionsparser.h"

#include <stdlib.h>

int main(int argc, char **argv){
	const char *filename = "src/microstran/properties.txt";
	if(argc>=2){
		filename = argv[1];
	}

	FILE *f;
	f = fopen(filename,"r");
	if(f==NULL){
		fprintf(stderr,"Unable to open section library '%s'",filename);
		exit(1);
	}

	const char *sn = "12ROD";
	if(argc >=3){
		sn = argv[2];
	}

	parse *p;
	p = parseCreateFileName(filename);

	section_library *l = NULL;
	l = section_library_create();
	parseSections(p,l);

	fprintf(stderr,"\n\nLooking up section '%s'...\n",sn);
	const section *s;
	s = section_find(l,sn);
	if(s==NULL){
		fprintf(stderr,"SECTION NOT FOUND\n");
		return 1;
	}else{
		fprintf(stderr,"SECTION FOUND\n");

		section_print(stderr,s);

		if(section_is_chs(s)){
			fprintf(stderr,"CHS\n");
			fprintf(stderr,"outside diameter = %f\n", section_chs_outside_diameter(s));
			fprintf(stderr,"thickness = %f\n\n", section_chs_thickness(s));
		}

		section_outline *o;
		o = section_get_outline(s);
		if(o){
			fprintf(stderr,"\nOutline:\n");
			section_outline_print(stderr,o);
		}
	}
	return 0;
}

