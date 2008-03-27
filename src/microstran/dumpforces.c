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
	little program to demonstrate that we can meaninfully parse model forces data files.
*/

#include "forceparser.h"

#include <stdlib.h>

int main(int argc, char **argv){
	const char *filename = "forces.p1";
	if(argc==2){
		filename = argv[1];
	}

	FILE *f;
	f = fopen(filename,"r");
	if(f==NULL){
		fprintf(stderr,"Unable to member forces data file '%s'",filename);
		exit(1);
	}

	parse *p;
	p = parseCreateFileName(filename);

	modelforces *mf = NULL;
	parseMicrostranForces(p,&mf);

	unsigned caseid = 1;
	fprintf(stderr,"Looking up case %d...\n",caseid);
	const caseforces *cf;
	cf = modelforces_find_case(mf,caseid);
	if(cf==NULL){
		fprintf(stderr,"LOAD CASE %d NOT FOUND\n",caseid);
		return 1;
	}else{
		fprintf(stderr,"CASEID %d FOUND\n",caseid);
	}
	return 0;
}

