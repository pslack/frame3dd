/** @FILE
	parser for microstran .arc files
*/

#define MSTRANP_BUILD

#include "displacementparser.h"

#include <string.h>

cbool parseNonHeading(parse *p){
	char c;
	char s[500];

	if(!parseAChar(p,&c)){
		return 0; /* end of document? */
	}

	/* line starting with '==' is a heading */

	if(c=='='){
		if(!parseAChar(p,&c)){
			parseUnParseChar(p,'=');
			//fprintf(stderr,"NONHEADING\n");
			return 0; /* failed to read two chars. unparse first '=' */
		}
		if(c=='='){
			/* we got two '==': unparse them and return fals */
			parseUnParseChar(p,'=');
			parseUnParseChar(p,'=');
			//fprintf(stderr,"NONHEADING\n");
			return 0;
		}
	}

	/* got something other than a heading. read to end of line */

	return (
		parseStrExcept(p,"\n\r",s,500)
		&& parseEOLplus(p)
		//&& assign(fprintf(stderr,"...%s\n",s))
	);
}
		
cbool parseHeading(parse *p, char *heading, unsigned maxlen){
	
	char c1, c2;
	char s[500];
	char *i, *j;
	if(!parseAChar(p,&c1)){
		return 0; /* end of document? */
	}

	if(c1!='='){
		parseUnParseChar(p,'=');
		return 0; /* not a '=' so not a heading */
	}

	if(!parseAChar(p,&c2)){
		parseUnParseChar(p,'=');
		return 0; /* failed to parse two chars */
	}

	if(c2!='='){
		/* we got two chars not both '=': unparse them and return fals */
		parseUnParseChar(p,c2);
		parseUnParseChar(p,c1);
		return 0;
	}

	//fprintf(stderr,"HEADING\n");
		
	if(!(maybe(parseWS(p))
		&& parseStrExcept(p,"\n\r",s,500)
		&& parseEOLplus(p)
		//&& assign(fprintf(stderr,"%s\n",s))
	)){
		fprintf(stderr,"FAILED HEADING MATCH\n");
		return 0;
	}

	for(i=s, j=heading; *i!='=' && *i!='\0' && j-heading < maxlen - 1; ++i, ++i, ++j){
		*j=*i;
	}
	*j='\0';

	//fprintf(stderr,"Heading '%s'\n",heading);
	if(j-heading+1==maxlen){
		fprintf(stderr,"TOO LONG\n");
		return 0;
	}
	return 1;
}

cbool parseToHeading(parse *p,const char *heading){
#define MAXHEADINGNAME 500
	char foundheading[MAXHEADINGNAME];
	cbool foundcorrect = 0;
	
	while(many(parseNonHeading(p)) && parseHeading(p,foundheading,MAXHEADINGNAME)){
		if(strcmp(foundheading,heading)==0){
			foundcorrect = 1;
			break;
		}
	}
	return foundcorrect;
}

cbool parseDisplacementNode(parse *p, casedisplacements *cd, const unsigned caseid){
	unsigned nodeid;
	double dx,dy,dz;
	double rx,ry,rz;
	return (
		parseNumber(p,&nodeid)
		//&& assign(fprintf(stderr,"NODE %d for case %d\n",nodeid,cd->id))
		&& parseWS(p)
		&& parseDouble(p,&dx)
		&& parseWS(p)
		&& parseDouble(p,&dy)
		&& parseWS(p)
		&& parseDouble(p,&dz)
		&& parseWS(p)
		&& parseDouble(p,&rx)
		&& parseWS(p)
		&& parseDouble(p,&ry)
		&& parseWS(p)
		&& parseDouble(p,&rz)
		&& parseEOLplus(p)
		&& casedisplacements_add_node(cd,nodeid,dx,dy,dz,rx,ry,rz)
	);
}

cbool parseDisplacementCASE(parse *p, modeldisplacements *md){
	unsigned caseid;
	char c;
	char casename[MAXCASENAME];
	casedisplacements *cd = NULL;
	return (
		parseThisString(p,"CASE")
		&& parseWS(p)
		&& parseNumber(p,&caseid) //&& assign(fprintf(stderr,"CASE %d",caseid))
		&& parseThisChar(p,':')
		&& parseWS(p)
		&& parseStrExcept(p,"\n\r",casename,MAXCASENAME) //&& assign(fprintf(stderr," '%s'\n",casename))
		&& parseEOLplus(p)
		&& parseThisString(p,"Node")
		&& many(parseCharExcept(p,"\n\r",&c))
		&& parseEOLplus(p)
		&& parseThisString(p,"m")
		&& parseWS(p)
		&& parseThisString(p,"m")
		&& parseWS(p)
		&& parseThisString(p,"m")
		&& parseWS(p)
		&& parseThisString(p,"rad")
		&& parseWS(p)
		&& parseThisString(p,"rad")
		&& parseWS(p)
		&& parseThisString(p,"rad")
		&& parseEOLplus(p)
		&& assign(cd = casedisplacements_create(caseid,casename))
		&& many(parseDisplacementNode(p, cd, caseid))
		//&& assign(fprintf(stderr,"Case %d contains %d nodes\n",caseid,ARRAY_NUM(cd->nodes)))
		&& modeldisplacements_add_case(md, cd)
	) || (
		(cd && assign(casedisplacements_destroy(cd)))
		&& fail
	);
}
		
cbool parseMicrostranDisplacements(parse *p, modeldisplacements **md){
	modeldisplacements *md1 = NULL;
	*md = NULL;
	return (
		assign(md1 = modeldisplacements_create())
		&& parseToHeading(p,"LOAD CASES")
		&& many(
			parseToHeading(p,"NODE DISPLACEMENTS")
			&& parseDisplacementCASE(p,md1)
		)
		&& assign(fprintf(stderr,"Displacement file contains %d cases\n",ARRAY_NUM(md1->cases)))
		&& assign(*md = md1)
		&& assign(fprintf(stderr,"Displacement file contains %d cases\n",ARRAY_NUM((*md)->cases)))
	) || (
		(
			md1 
			&& (modeldisplacements_destroy(md1),1)
		)
		&& fail
	);
	
	// read file header
	// skip down to 'N O D E   D I S P'...
	// read in list of node IDs and displacements
	// move to next load case
}

#if 0
model *applyMicrostranDisplacements(model *m, modeldisplacements **d){
	// create a copy of the mode with the displacements applied.
	// check that displacements are given for ALL nodes in the model
	// output the maximum displacement magnitude ccc
}
#endif


