/** @FILE
	parser for our home-made sections database files (eg 'properties.txt')
*/

#define MSTRANP_BUILD

#include <stdio.h>
#include <string.h>
#include "parse.h"
#include "error.h"
#include "case.h"
#include "sectionsparser.h"

cbool parseSectionsFileComment(parse *p){
	char c;
	return (
		parseThisChar(p, '#')/* && assign(fprintf(stderr,"#")) */
		&& many(parseCharExcept(p,"\n\r",&c)/* && fputc(c,stderr)*/)
		&& parseEOLplus(p) /*&& assign(fprintf(stderr,"\n"))*/
	);
}

/**
	Parse the fields for a CHS section (we've already read the name)
*/
cbool parseCHS(parse *p, section *s1){
	struct section_chs_struct s;
	double ignore;
	return (
		maybe(parseWS(p))
		&& parseThisString(p,"CHS")/* && assign(fprintf(stderr,"CHS"))*/
		&& parseWS(p)
		&& parseDouble(p,&s.d_o)// && assign(fprintf(stderr,"\td_o = %f\n",s.d_o))
		&& parseWS(p)
		&& parseDouble(p, &s.t)// && assign(fprintf(stderr,"\tt = %f\n",s.t))
		&& parseWS(p)
		&& parseDouble(p, &ignore)// && assign(fprintf(stderr,"\tignore %lf\n",ignore))
		&& parseWS(p)
		&& parseDouble(p, &ignore)// && assign(fprintf(stderr,"\tignore %lf\n",ignore))
		&& parseWS(p)
		&& parseDouble(p, &ignore)// && assign(fprintf(stderr,"\tignore %lf\n",ignore))
		&& parseWS(p)
		&& parseDouble(p, &ignore)// && assign(fprintf(stderr,"\tignore %lf\n",ignore))
		&& parseWS(p)
		&& parseDouble(p,&s.A)// && assign(fprintf(stderr,"\tA = %lf\n",s.A))
		&& parseWS(p)
		&& parseDouble(p,&s.I)// && assign(fprintf(stderr,"\tI = %lf\n",s.I))
		&& parseWS(p)
		&& parseDouble(p,&s.Z)// && assign(fprintf(stderr,"\tZ = %lf\n",s.Z))
		&& parseWS(p)
		&& parseDouble(p,&s.S)// && assign(fprintf(stderr,"\tS = %lf\n",s.S))
		&& parseWS(p)
		&& parseDouble(p,&s.r)
		&& parseWS(p)
		&& parseDouble(p,&s.J)
		&& parseWS(p)
		&& parseDouble(p,&s.C)
		&& parseWS(p)
		&& parseDouble(p,&s.k_f)
		&& parseWS(p)
		&& parseDouble(p,&s.lambda_s)// && assign(fprintf(stderr,"\tlambda_s = %lf\n",s.lambda_s))
		&& parseWS(p)
		&& parseAChar(p,&s.compactness)// && assign(fprintf(stderr,"\tcompactness = %c\n",s.compactness))
		&& parseWS(p)
		&& parseDouble(p,&s.Ze)/* && assign(fprintf(stderr,"\tZe = %lf\n",s.Ze))*/
		&& maybe(parseWS(p))
		&& parseEOLplus(p)
		&& assign(s1->type = SECTION_CHS)
		&& assign(s1->chs = s)
	);
}

/**
	Parse the fields for a I-section (we've already read the name)
*/
cbool parseIsec(parse *p, section *s1){
	struct section_isec_struct s;
	double ignore;
	return (
		maybe(parseWS(p))
		&& ((parseThisString(p,"UB") /* && assign(fprintf(stderr,"UB"))*/)
			|| (parseThisString(p,"UC") /* && assign(fprintf(stderr,"UC"))*/))
		&& parseWS(p)
		&& parseDouble(p,&s.m_linear)
		&& parseWS(p)
		&& parseDouble(p, &s.d)
		&& parseWS(p)
		&& parseDouble(p, &s.bf)
		&& parseWS(p)
		&& parseDouble(p, &s.tf)
		&& parseWS(p)
		&& parseDouble(p, &s.tw)
		&& parseWS(p)
		&& parseDouble(p, &s.r1)
		&& parseWS(p)
		&& parseDouble(p,&s.Ag)
		&& parseWS(p)
		&& parseDouble(p,&s.lx)
		&& parseWS(p)
		&& parseDouble(p,&s.Zx)
		&& parseWS(p)
		&& parseDouble(p,&s.Sx)
		&& parseWS(p)
		&& parseDouble(p,&s.rx)
		&& parseWS(p)
		&& parseDouble(p,&s.ly)
		&& parseWS(p)
		&& parseDouble(p,&s.Zy)
		&& parseWS(p)
		&& parseDouble(p,&s.Sy)
		&& parseWS(p)
		&& parseDouble(p,&s.ry)
		&& parseWS(p)
		&& parseDouble(p,&s.J)
		&& parseWS(p)
		&& parseDouble(p,&s.lw)
		&& maybe(parseWS(p))
		&& parseEOLplus(p)
		&& assign(s.d *= 1e-3)
		&& assign(s.bf *= 1e-3)
		&& assign(s.tf *= 1e-3)
		&& assign(s.tw *= 1e-3)
		&& assign(s.r1 *= 1e-3)
		&& assign(s.Ag *= 1e-6)

		&& assign(s1->type = SECTION_ISEC)
		&& assign(s1->isec = s)
	);
}

cbool parseSection(parse *p, section_library *l){
	section s;
	return (
		parseStrExcept(p," \n\r\t",s.name,SECTION_NAME_MAX) && strlen(s.name)// && assign(fprintf(stderr,"Section \"%s\"",s.name))
		&& (
			parseCHS(p, &s)
			|| parseIsec(p, &s)
			/* || parseUB(...) */
			/* || parseRHS(...) */
			/* etc */
		) && assign(array_append(&(l->a),&s)) 
		//&& assign(section_print(stderr,&s))
	);
}

/* SECTIONS LIBRARY FILE*/

cbool parseSections(parse *p, section_library *l){
	return (
		many(
			parseSectionsFileComment(p)
			|| parseSection(p,l)
		)
	);
}


