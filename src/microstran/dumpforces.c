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

	modelforces *MF = NULL;
	parseMicrostranForces(p,&MF);
	if(MF==NULL){
		fprintf(stderr,"ERROR: Failed to parse forces!\n");
		exit(1);
	}

	unsigned caseid = 1;
	fprintf(stderr,"Looking up case %d...\n",caseid);
	caseforces *cf;
	cf = modelforces_find_case(MF,caseid);
	if(cf==NULL){
		fprintf(stderr,"LOAD CASE %d NOT FOUND\n",caseid);
		return 1;
	}else{
		fprintf(stderr,"CASE ID %d FOUND\n",caseid);
	}

	unsigned memberid = 4;
	const memberforce *mf;
	mf = caseforces_find_member(cf, memberid);

	if(mf==NULL){
		fprintf(stderr,"MEMBER FORCE FOR MEMBER %d NOT FOUND\n",memberid);
		return 1;
	}
	fprintf(stderr,"FOUND FORCE ON MEMBER %d\n",memberid);
	
	fprintf(stderr,"Force at 'from' end (node %d):\n",mf->fromnode.node);
	VEC3_PR(mf->fromnode.F);
	VEC3_PR(mf->fromnode.M);

	fprintf(stderr,"Force at 'to' end (node %d):\n",mf->tonode.node);
	VEC3_PR(mf->tonode.F);
	VEC3_PR(mf->tonode.M);
	
	return 0;
}

