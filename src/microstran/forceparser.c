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
*//** @FILE
	parser for microstran .p1 member forces data files

	@NOTE there are important differences between the axes and sign conventions
	used here and those specified by the Microstran .arc format. @ENDNOTE
*/

#define MSTRANP_BUILD

#include "forceparser.h"

#include <string.h>

static cbool parseNonHeading(parse *p){
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

static cbool parseHeading(parse *p, char *heading, unsigned maxlen){

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

static cbool parseToHeading(parse *p,const char *heading){
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

cbool parseMemberNodeForce(parse *p, membernodeforce *mnf){
	unsigned nodeid;
	vec3 F, M;
	return (
		//assign(fprintf(stderr,"Looking for member node force...\n")) &&
		parseNumber(p,&nodeid)
		//&& assign(fprintf(stderr,"NODE %d\n",nodeid))
		&& parseWS(p)
		&& parseDouble(p,&F.z)
		&& parseWS(p)
		&& parseDouble(p,&F.y)
		&& parseWS(p)
		&& parseDouble(p,&F.x)
		//&& assign(VEC3_PR(F))
		&& parseWS(p)
		&& parseDouble(p,&M.z)
		&& parseWS(p)
		&& parseDouble(p,&M.y)
		&& parseWS(p)
		&& parseDouble(p,&M.x)
		//&& assign(VEC3_PR(M))
		&& parseEOLplus(p)
		&& assign(F.x = -F.x)
		&& assign(M.z = -M.z)
		&& assign(*mnf = membernodeforce_create(nodeid,F,M))
		//&& assign(fprintf(stderr,"Member node force on node %d OK\n",nodeid))
	) || (
		assign(fprintf(stderr,"Failed to parse member node force\n"))
		&& fail
	);
}

cbool parseMemberForce(parse *p, caseforces *cf, const unsigned caseid){
	unsigned memberid = 0;
	membernodeforce fromnode;
	membernodeforce tonode;
	return (
		parseNumber(p,&memberid)
		//&& assign(fprintf(stderr,"Member %d...\n",memberid))
		&& parseWS(p)
		&& parseMemberNodeForce(p, &fromnode)
		//&& assign(fprintf(stderr,"Got 'from' node for member %d...\n",memberid))
		//&& assign(fprintf(stderr,"Looking for WS...\n"))
		&& parseMemberNodeForce(p, &tonode)
		//&& assign(fprintf(stderr,"Got 'to' node for member %d...\n",memberid))
		&& caseforces_add_member(cf, memberid, fromnode, tonode)
		//&& assign(fprintf(stderr,"Added member %d to case %d OK...\n",memberid,caseid))
		//&& assign(fprintf(stderr,"Member %d force OK\n",memberid))
	)/* || (
		assign(fprintf(stderr,"Failed to parse forces for member %d\n",memberid))
		&& fail
	)*/;
}

cbool parseForceCASE(parse *p, modelforces *mf){
	unsigned caseid;
	char c;
	char casename[MAXCASENAME];
	caseforces *cf = NULL;
	return (
		parseThisString(p,"CASE")
		&& parseWS(p)
		&& parseNumber(p,&caseid) //&& assign(fprintf(stderr,"CASE %d",caseid))
		&& parseThisChar(p,':')
		&& parseWS(p)
		&& parseStrExcept(p,"\n\r",casename,MAXCASENAME) //&& assign(fprintf(stderr," '%s'\n",casename))
		&& parseEOLplus(p)
		&& parseThisString(p,"Member")
		&& many(parseCharExcept(p,"\n\r",&c))
		&& parseEOLplus(p)
		&& parseThisString(p,"kN")
		&& parseWS(p)
		&& parseThisString(p,"kN")
		&& parseWS(p)
		&& parseThisString(p,"kN")
		&& parseWS(p)
		&& parseThisString(p,"kNm")
		&& parseWS(p)
		&& parseThisString(p,"kNm")
		&& parseWS(p)
		&& parseThisString(p,"kNm")
		&& parseEOLplus(p)
		&& assign(cf = caseforces_create(caseid,casename))
		&& one_or_more(parseMemberForce(p, cf, caseid))
		//&& assign(fprintf(stderr,"Case %d contains %d nodes\n",caseid,ARRAY_NUM(cf->nodes)))
		&& modelforces_add_case(mf, cf)
	) || (
		(cf && assign(caseforces_destroy(cf)))
		&& fail
	);
}

cbool parseMicrostranForces(parse *p, modelforces **mf){
	modelforces *mf1 = NULL;
	*mf = NULL;
	return (
		assign(mf1 = modelforces_create())
		&& parseToHeading(p,"LOAD CASES")
		&& one_or_more(
			parseToHeading(p,"MEMBER FORCES")
			&& parseForceCASE(p,mf1)
		)
		&& assign(*mf = mf1)
		&& assign(fprintf(stderr,"Member forces data file parsed; contains %d load-cases.\n",ARRAY_NUM((*mf)->cases)))
	) || (
		(
			mf1
			&& (modelforces_destroy(mf1),1)
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


