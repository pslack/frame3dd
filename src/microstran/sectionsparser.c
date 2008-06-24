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
	parser for our home-made sections database files (eg 'properties.txt')
*/

#define MSTRANP_BUILD

#include <stdio.h>
#include <string.h>
#include "parse.h"
#include "error.h"
#include "case.h"
#include "sectionsparser.h"

#define PI	3.141592653589793

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

/*
#Dimesions		-------------------------------------------------------------------------------PROPERTIES-------------------------------------------------------------------------------------------------------------------											-------- PROPERTIES FOR DESIGN TO AS 4100 -------					
#Designation		Mass	External		Gross		About x-, y- and n-axis				Torsion	Torsion	Form	About x- and y-axis 				
#		per m	Surface Area		section area		----------------------------------------------------------------------				Constant	Modulus	Factor	----------------------------------------------------				
#d x b x t		m_linear	per_m	per_t	Ag	Ix	Zx	Zn	Sx	rx	J	C	kf	lambda_e	Compactness	Ze		
#mm x mm x mm		kg/m	m²/m	m²/t	mm²	10^6mm^4	10³mm³	10³mm³	10³mm³	mm	10^6mm^4	10³mm³			(C,N,S)	10³mm³		
250x250x9.0SHS	SHS	65.9	0.96	14.6	8400	79.8	639	477	750	97.5	129	972	1	30.5	N	744		
250x250x6.0SHS	SHS	45	0.97	21.7	5730	56.2	450	330	521	99	88.7	681	0.85	46.9	S	409		
*/
cbool parseSHS(parse *p, section *s1){
	struct section_shs_struct s;
	double ignore;
	return (
#define ASSIGN(Z) 1
		maybe(parseWS(p))
		&& parseThisString(p,"SHS") && ASSIGN(fprintf(stderr,"SHS"))
		&& parseWS(p)
		&& parseDouble(p,&s.d) && ASSIGN(fprintf(stderr,"\td = %f\n",s.d))
		&& parseWS(p)
		&& parseDouble(p, &s.t) && ASSIGN(fprintf(stderr,"\tt = %f\n",s.t))
		&& parseWS(p)
		&& parseDouble(p, &s.m_linear) && ASSIGN(fprintf(stderr,"\tm_linear = %lf\n",s.m_linear))
		&& parseWS(p)
		&& parseDouble(p, &ignore)// && assign(fprintf(stderr,"\tignore %lf\n",ignore))
		&& parseWS(p)
		&& parseDouble(p, &ignore)// && assign(fprintf(stderr,"\tignore %lf\n",ignore))
		&& parseWS(p)
		&& parseDouble(p, &s.Ag) && ASSIGN(fprintf(stderr,"\tAg = %lf\n",s.Ag))
		&& parseWS(p)
		&& parseDouble(p,&s.Ix) && ASSIGN(fprintf(stderr,"\tIx = %lf\n",s.Ix))
		&& parseWS(p)
		&& parseDouble(p,&s.Zx) && ASSIGN(fprintf(stderr,"\tZx = %lf\n",s.Zx))
		&& parseWS(p)
		&& parseDouble(p,&s.Zn) && ASSIGN(fprintf(stderr,"\tZn = %lf\n",s.Zn))
		&& parseWS(p)
		&& parseDouble(p,&s.Sx) && ASSIGN(fprintf(stderr,"\tSx = %lf\n",s.Sx))
		&& parseWS(p)
		&& parseDouble(p,&s.rx) && ASSIGN(fprintf(stderr,"\trx = %lf\n",s.rx))
		&& parseWS(p)
		&& parseDouble(p,&s.J) && ASSIGN(fprintf(stderr,"\tJ = %lf\n",s.J))
		&& parseWS(p)
		&& parseDouble(p,&s.C) && ASSIGN(fprintf(stderr,"\tc = %lf\n",s.C))
		&& parseWS(p)
		&& parseDouble(p,&s.k_f) && ASSIGN(fprintf(stderr,"\tk_f = %lf\n",s.k_f))
		&& parseWS(p)
		&& parseDouble(p,&s.lambda_e)// && assign(fprintf(stderr,"\tlambda_s = %lf\n",s.lambda_s))
		&& parseWS(p)
		&& parseAChar(p,&s.compactness)// && assign(fprintf(stderr,"\tcompactness = %c\n",s.compactness))
		&& parseWS(p)
		&& parseDouble(p,&s.Ze)/* && assign(fprintf(stderr,"\tZe = %lf\n",s.Ze))*/
		&& maybe(parseWS(p))
		&& parseEOLplus(p)
		&& assign(s.d *= 1e-3)
		&& assign(s.t *= 1e-3)
		&& assign(s.rx *= 1e-3)
		&& assign(s.Ag *= 1e-6)
		&& assign(s1->type = SECTION_SHS)
		&& assign(s1->shs = s)
#undef ASSIGN
	);
}


/**
	Parse the fields for a Top Hat section (we've already read the name)
*/
cbool parseTopHat(parse *p, section *s1){
	struct section_tophat_struct s;
	s.inverted = 0;
	return (
		maybe(parseWS(p))
		&& (parseThisString(p,"TOPHAT"))
		&& parseWS(p)
		&& parseDouble(p,&s.a)
		&& parseWS(p)
		&& parseDouble(p, &s.b)
		&& parseWS(p)
		&& parseDouble(p, &s.c)
		&& parseWS(p)
		&& parseDouble(p, &s.d)
		&& parseWS(p)
		&& parseDouble(p, &s.theta)
		&& parseWS(p)
		&& parseDouble(p, &s.phi)
		&& parseWS(p)
		&& parseDouble(p,&s.t)
		&& maybe(
			parseWS(p)
			&& (
				parseThisString(p,"INVERTED")
				&& assign(s.inverted = 1)
			)
		) && parseEOLplus(p)
		&& assign(s.a *= 1e-3)
		&& assign(s.b *= 1e-3)
		&& assign(s.c *= 1e-3)
		&& assign(s.d *= 1e-3)
		&& assign(s.theta *= PI/180.)
		&& assign(s.phi *= PI/180.)
		&& assign(s.t *= 1e-3)
		&& assign(s1->type = SECTION_TOPHAT)
		&& assign(s1->tophat = s)
	);
}

/**
	Parse the fields for a rod section (we've already read the name)
*/
cbool parseRod(parse *p, section *s1){
	struct section_rod_struct s;
	return (
		maybe(parseWS(p))
		&& parseThisString(p,"ROD")
		&& parseWS(p)
		&& parseDouble(p,&s.d)
		&& parseWS(p)
		&& parseDouble(p,&s.A)
		&& parseWS(p)
		&& parseDouble(p,&s.Ax)
		&& parseWS(p)
		&& parseDouble(p,&s.J)
		&& parseWS(p)
		&& parseDouble(p,&s.Ix)
		&& parseWS(p)
		&& parseDouble(p,&s.M)
		&& parseEOLplus(p)
		&& assign(s.d *= 1e-3)
		&& assign(s1->type = SECTION_ROD)
		&& assign(s1->rod = s)
	);
}	


cbool parseSection(parse *p, section_library *l){
	section s;
#define ASSIGN(Z) 1
	return (
		parseStrExcept(p," \n\r\t",s.name,SECTION_NAME_MAX) && strlen(s.name)// && assign(fprintf(stderr,"Section \"%s\"",s.name))
		&& (
			parseCHS(p, &s)
			|| parseIsec(p, &s)
			|| parseSHS(p, &s)
			|| parseTopHat(p, &s)
			|| parseRod(p, &s)
			/* || parseUB(...) */
			/* etc */
		) && assign(array_append(&(l->a),&s)) 
		&& ASSIGN(section_print(stderr,&s))
	);
#undef ASSIGN
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


