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
	little program to demonstrate that we can retrieve load case data from
	parsed microstran .arc files
*/

#include "modelparser.h"
#include "case.h"

#include <assert.h>
#include <stdlib.h>

int main(int argc, char **argv){
	const char *filename = "model.arc";
	if(argc==2){
		filename = argv[1];
	}

	FILE *f;
	f = fopen(filename,"r");
	if(f==NULL){
		fprintf(stderr,"Unable to .arc model file '%s'",filename);
		exit(1);
	}

	parse *p;
	p = parseCreateFileName(filename);

	model *M = NULL;
	parseModelMicrostran(p,&M);
	if(M==NULL){
		fprintf(stderr,"ERROR: Failed to parse model '%s'!\n",filename);
		exit(1);
	}

	fprintf(stderr,"======================================\n");

	unsigned i, j, nc = ARRAY_NUM(M->cases);
	for(i=0; i<nc; ++i){
		case_stmt *c = array_get(&M->cases, i);
		switch(c->type){
			case CASE_LOADS:
				fprintf(stderr,"Case %d: %d NDLDs '%s'\n", c->id, ARRAY_NUM(c->data), c->name);
				break;
			case CASE_COMB:
				fprintf(stderr,"\nCase %d: %d COMBs '%s'\n", c->id, ARRAY_NUM(c->data), c->name);
				for(j = 0; j<ARRAY_NUM(c->data); ++j){
					comb_stmt *c1 = array_get(&c->data, j);
					case_stmt *case1 = array_get(&M->cases, c1->caseindex);
					//fprintf(stderr,"c1->c = %p, M->cases.data = %p\n", case1, M->cases.data);
					//assert(c1->c >= (case_stmt *)M->cases.data);
					//assert(c1->c < (((case_stmt *)M->cases.data + ARRAY_NUM(M->cases))));
					fprintf(stderr," + %f × CASE %d '%s'\n", c1->factor, case1->id, case1->name);
				}
				break;
			default:
				fprintf(stderr,"CASE UNKNOWN TYPE\n");
				exit(1);
		}
	};
	
	return 0;
}

