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
#include <stdio.h>
#include "modelparser.h"
#include "displacementparser.h"
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
#if 0
	fprintf(stderr,"STARTING PARSE...\n");
	parse *p;
	model *a; /* will be allocated by the parser */
	p = parseCreateFile(stdin);
	if(parseMicrostranArchive(p,&a)){
		fprintf(stderr,"\nParsed OK\n%d nodes, %d members, %d section properties, %d materials read.\n"
			,a->num_nodes,a->num_membs,a->num_props,a->num_matls
		);
	}else{
		fprintf(stderr,"\n\nPARSE FAILED\n");
	}
	return 0;
#else
	parse *p;

	FILE *modelfile;
	const char *defaultmodelfname = "model.arc";

	FILE *displacementfile;
	const char *defaultdisplacementfname = "model.p1";

	unsigned caseid;
	const unsigned defaultcaseid = 1;

	const char *modelfname;
	const char *displacementfname;

	if(argc>1){
		modelfname = argv[1];
		fprintf(stderr,"Model file '%s'\n",modelfname);
	}else{
		modelfname = defaultmodelfname;
		fprintf(stderr,"Using default model filename '%s'.\n",modelfname);
	}

	modelfile = fopen(modelfname,"r");

	if(!modelfile){
		fprintf(stderr,"\nFailed to read model file '%s'\n",modelfname);
		return 1;
	}else{
		fprintf(stderr,"Opened model file '%s'...\n",modelfname);
	}

	if(argc>2){
		displacementfname = argv[2];
		fprintf(stderr,"Displacement file '%s'\n",displacementfname);
	}else{
		displacementfname = defaultdisplacementfname;
		fprintf(stderr,"Using default displacement filename '%s'.\n",displacementfname);
	}

	displacementfile = fopen(displacementfname,"r");
	if(!displacementfile){
		fprintf(stderr,"\nFAILED TO OPEN DISPLACEMENT FILE '%s'\n",displacementfname);
		return 1;
	}else{
		fprintf(stderr,"Opened displacement file '%s'...\n",displacementfname);
	}

	fprintf(stderr,"Loading model...\n");
	p = parseCreateFile(modelfile);
	if(!p){
		fprintf(stderr,"\nERROR: Failed to load model file '%s'\n", modelfile);
		return 1;
	}
	model *m;
	if(parseModelMicrostran(p,&m)){
		fprintf(stderr,"\nParsed OK\n%d nodes, %d members, %d section properties, %d materials read.\n"
			,m->num_nodes,m->num_membs,m->num_props,m->num_matls
		);
	}else{
		fprintf(stderr,"\n\nPARSE FAILED\n");
		return 1;
	}
	parseDispose(p);

	fprintf(stderr,"Loading displacements...\n");
	parse *p1;
	p1 = parseCreateFile(displacementfile);
	modeldisplacements *md;
	if(parseMicrostranDisplacements(p1,&md)){
		fprintf(stderr,"Parse OK\n");
	}else{
		fprintf(stderr,"\n\nPARSE FAILED\n");
	}

	if(argc>3){
		caseid = atoi(argv[3]);
		fprintf(stderr,"Applying displacements for load case %d...\n",caseid);
	}else{
		caseid = defaultcaseid;
		fprintf(stderr,"Applying default load case %d.\n",caseid);
	}

	casedisplacements *cd;
	cd = modeldisplacements_find_case(md, caseid);
	if(!cd){
		fprintf(stderr,"\n\nFAILED TO FIND CASE %d\n",caseid);
		return 1;
	}

	if(!model_apply_displacements(m,cd)){
		fprintf(stderr,"\n\nFAILED TO APPLY DISPLACEMENTS\n");
	}
	return 0;
#endif
}
