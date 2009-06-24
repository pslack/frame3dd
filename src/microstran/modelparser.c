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
	parser for microstran .arc files
*/

#define MSTRANP_BUILD

#include <stdio.h>
#include <string.h>
#include "parse.h"
#include "error.h"
#include "case.h"
#include "modelparser.h"

cbool parseComments(parse *p){
	char c;
	return many(
		parseThisChar(p, '*') && assign(fprintf(stderr,"*"))
		&& many(parseCharExcept(p,"\n\r",&c) && fputc(c,stderr))
		&& parseEOLplus(p) && assign(fprintf(stderr,"\n"))

	);
}

/*----------------------------------------------------------------------------*/

/* line with keyword followed by a number */
cbool parseNumericLine(parse *p, const char *kw, unsigned *valptr){
	return (
		parseComments(p)
		&& parseThisString(p,kw)
		&& parseWS(p)
		&& parseNumber(p, valptr)
		&& maybe(parseWS(p))
		&& parseEOLplus(p)
	);
}


/* UNIT statement */

cbool parseUNIT(parse *p, unit_stmt *u){
	return (
		parseComments(p)
		&& parseThisString(p,"UNIT")
		&& parseWS(p)
		&& parseNumber(p,&(u->num))
		&& parseWS(p)
		&& parseNonWS(p,u->lengthunit,MAXUNIT)
			//&& assign(fprintf(stderr,"LENGTH UNIT %s\n",u->lengthunit))
		&& parseWS(p)
		&& parseNonWS(p,u->forceunit,MAXUNIT)
			//&& assign(fprintf(stderr,"FORCE UNIT %s\n",u->forceunit))
		&& parseWS(p)
		&& parseNonWS(p,u->massunit,MAXUNIT)
			//&& assign(fprintf(stderr,"MASS UNIT %s\n",u->massunit))
		&& parseWS(p)
		&& parseNonWS(p,u->tempunit,MAXUNIT)
			//&& assign(fprintf(stderr,"TEMP UNIT %s\n",u->tempunit))
		&& parseEOLplus(p)
	);
}

/* NODE statement */

cbool parseNODE(parse *p, model *a){
	unsigned nodeid, flags = 0;
	vec3 pos;
	return (
		parseComments(p)
		&& parseThisString(p,"NODE")
		&& parseWS(p)
		&& parseNumber(p,&nodeid)
		&& parseWS(p)
		&& parseDouble(p, &(pos.x))
		&& parseWS(p)
		&& parseDouble(p, &(pos.y))
		&& parseWS(p)
		&& parseDouble(p, &(pos.z))
		&& parseWS(p)
		&& parseBitChar(p,MSTRANP_NODE_FIXX, &flags)
		&& parseBitChar(p,MSTRANP_NODE_FIXY, &flags)
		&& parseBitChar(p,MSTRANP_NODE_FIXZ, &flags)
		&& parseBitChar(p,MSTRANP_NODE_FIXMX, &flags)
		&& parseBitChar(p,MSTRANP_NODE_FIXMY, &flags)
		&& parseBitChar(p,MSTRANP_NODE_FIXMZ, &flags)
		&& parseEOLplus(p)
		&& model_add_node(a,node_create(nodeid,pos,flags))
		//&& assign(fprintf(stderr,"NODE %d\n",nodeid))
	);
}

/* MEMB statement */

cbool parseMemberOrientation(parse *p, model *a, member_orientation *orient){
	orient->axis = '\0';
	orient->dir = '+';
	return(
		(
			maybe(parseCharIn(p,"-+",&(orient->dir)) /*assign(fprintf(stderr,"%c",*axisdir))*/)
			&& (
				(parseCharIn(p,"XYZ",&(orient->axis)) /*&& assign(fprintf(stderr,"%c",*axis))*/)
			)
		) || (
			//assign(fprintf(stderr,"aligntomemb="))
			parseNumber(p,&(orient->node)) //&& assign(fprintf(stderr,"%d",*aligntomemb))
			&& assign(orient->axis = '\0')
		)
	);
}

cbool parseMEMB(parse *p, model *a){
	unsigned id = 0;
	unsigned fromnode,tonode, prop, matl;
	member_orientation orient;
	unsigned flags1, flags2;
	return (
		parseComments(p)
		&& parseThisString(p,"MEMB") //&& assign(fprintf(stderr,"MEMB"))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&id) //&& assign(fprintf(stderr,"%d",id))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&fromnode) //&& assign(fprintf(stderr,"from=%d",fromnode))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&tonode) //&& assign(fprintf(stderr,"tonode=%d",tonode))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseMemberOrientation(p, a, &orient)
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&prop) //&& assign(fprintf(stderr,"prop=%d",prop))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&matl) //&& assign(fprintf(stderr,"matl=%d",matl))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&flags1) //&& assign(fprintf(stderr,"flags=%d",flags1))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&flags2) //&& assign(fprintf(stderr,"%d",flags2))
		&& parseEOLplus(p) //&& assign(fprintf(stderr,"\n"))
		&& model_add_memb(a,id,fromnode,tonode,orient,prop,matl,flags1,flags2)
		//&& assign(fprintf(stderr,"MEMB %d OK\n",id))
	);
}

/* MOFF statements: member offsets? */

cbool parseMOFF(parse *p, model *a){
	unsigned memberid;
	vec3 deltaA, deltaB;
	coord_sys_t coordsys;
	static int warnMOFF = 0;
	return (
		parseComments(p)
		&& parseThisString(p,"MOFF")
		&& parseWS(p)
		&& parseNumber(p,&memberid)
		&& parseWS(p)
		&& (
			(parseThisString(p,"LO") && assign(coordsys=MSTRANP_COORDS_LOCAL))
			|| (parseThisString(p,"GL") && assign(coordsys=MSTRANP_COORDS_GLOBAL))
		)
		&& parseWS(p)
		&& parseDouble(p,&deltaA.x)
		&& parseWS(p)
		&& parseDouble(p,&deltaA.y)
		&& parseWS(p)
		&& parseDouble(p,&deltaA.z)
		&& parseWS(p)
		&& parseDouble(p,&deltaB.x)
		&& parseWS(p)
		&& parseDouble(p,&deltaB.y)
		&& parseWS(p)
		&& parseDouble(p,&deltaB.z)
		&& parseEOLplus(p)
		&& model_add_member_offset(a, memberid, coordsys, deltaA, deltaB)
		&& assign(warnMOFF ? 0 : (fprintf(stderr,"WARNING: 'MOFF' implementation is still experimental.\n"), warnMOFF=1) )
	);
}

/* MSPR statements: member springs? */

// MSPR  8124  1500.000         R         R         R         R         R
cbool parseMSPR(parse *p, model *a){
	unsigned memberid;
	double stiffness;
	static int warnMSPR = 0;
	return (
		parseComments(p)
		&& parseThisString(p,"MSPR")
		&& parseWS(p)
		&& parseNumber(p,&memberid) //&& assign(fprintf(stderr,"MSPR %d\n",memberid))
		&& parseWS(p)
		&& parseDouble(p,&stiffness)
		&& parseWS(p)
		&& parseThisString(p,"R")
		&& parseWS(p)
		&& parseThisString(p,"R")
		&& parseWS(p)
		&& parseThisString(p,"R")
		&& parseWS(p)
		&& parseThisString(p,"R")
		&& parseWS(p)
		&& parseThisString(p,"R")
		&& parseEOLplus(p)
		//&& assign(fprintf(stderr,"MSPR %d\n",memberid))
		&& assign(warnMSPR ? 0 : (fprintf(stderr,"WARNING: 'MSPR' statements ignored.\n"), warnMSPR=1) )

	);
}

/* MTYP statements: member types (eg tension-only members) */
cbool parseMTYP(parse *p, model *a){
	// MTYP  6046  TONLY
	static int warnMTYP = 0;
	unsigned memberid;
	return (
		parseComments(p)
		&& parseThisString(p,"MTYP")
		&& parseWS(p)
		&& parseNumber(p,&memberid)
		&& parseWS(p)
		&& (
			parseThisString(p,"TONLY")
			|| parseThisString(p,"CONLY")
		)
		&& parseEOLplus(p)
		&& assign(warnMTYP ? 0 : (fprintf(stderr,"WARNING: 'MTYP' statements ignored.\n"),warnMTYP=1))
	);
}

/* PROP statement */

/**
	FIXME this should also parse LIBR lines with 'X' instead of 'Y' as the
	major-axis alignment
*/
cbool parsePROP(parse *p, model *a){
	unsigned id;
	char libr[MAXPROPLIBNAME] = "";
	char name[MAXPROPNAME] = "";
	char desc[MAXPROPDESC] = "";
	double vals[6];
#define ASSIGN(Z) 1
	return (
		parseComments(p)
		&& parseThisString(p,"PROP") && ASSIGN(fprintf(stderr,"PROP"))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseNumber(p,&id) && ASSIGN(fprintf(stderr,"%d",id))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& ((
			parseThisString(p,"LIBR") && ASSIGN(fprintf(stderr,"LIBR"))
			&& parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseStrExcept(p," \n\r\t",libr,MAXPROPLIBNAME) && ASSIGN(fprintf(stderr,"%s",libr))
			&& parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseStrExcept(p," \n\r\t",name,MAXPROPNAME) && ASSIGN(fprintf(stderr,"%s",name))
			&& parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseThisChar(p,'Y') && ASSIGN(fprintf(stderr,"Y"))
			&& maybe(
				parseWS(p)
				&& parseStrExcept(p,"\n\r\t",desc,MAXPROPDESC)
				&& ASSIGN(fprintf(stderr," %s",desc))
			) && parseEOLplus(p)
		) || (
			parseThisString(p,"PRIS") && ASSIGN(fprintf(stderr,"PRIS"))
			&& parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseStrExcept(p," \n\r\t",name,MAXPROPNAME) && ASSIGN(fprintf(stderr,"%s",name))
			&& parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseStrExcept(p,"\n\r\t",desc,MAXPROPDESC) && ASSIGN(fprintf(stderr,"%s",desc))
			&& parseEOLplus(p)
		))
		&& ASSIGN(fprintf(stderr,"\n\t"))
		&& parseDouble(p,vals+0) && ASSIGN(fprintf(stderr,"%e",vals[0]))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,vals+1) && ASSIGN(fprintf(stderr,"%e",vals[1]))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,vals+2) && ASSIGN(fprintf(stderr,"%e",vals[2]))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,vals+3) && ASSIGN(fprintf(stderr,"%e",vals[3]))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,vals+4) && ASSIGN(fprintf(stderr,"%e",vals[4]))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,vals+5) && ASSIGN(fprintf(stderr,"%e",vals[5]))
		&& parseEOLplus(p) && ASSIGN(fprintf(stderr,"\n"))
		&& model_add_prop(a, id, libr, name, desc, vals)
	);
#undef ASSIGN
}

/* MATL statement */

cbool parseMATL(parse *p, model *a){
	unsigned id;
	double E; /* force unit / length_unit^2 */
	double sigma_y; /* force unit / length_unit^2 */
	double rho; /* mass_unit / length_unit^3 */
	double beta; /* temp_unit^-1 ...guessing about this one */
#define ASSIGN(Z) 1
	return (
		parseComments(p)
		&& parseThisString(p,"MATL") && ASSIGN(fprintf(stderr,"MATL"))
		&& parseWS(p)
		&& parseNumber(p,&id) && ASSIGN(fprintf(stderr,"%d",id))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,&E) && ASSIGN(fprintf(stderr,"%e",E))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,&sigma_y) && ASSIGN(fprintf(stderr,"%e",sigma_y))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,&rho) && ASSIGN(fprintf(stderr,"%e",rho))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseDouble(p,&beta) && ASSIGN(fprintf(stderr,"%e",beta))
		&& parseEOLplus(p) && ASSIGN(fprintf(stderr,"\n"))
		&& model_add_matl(a,id,E,sigma_y,rho,beta)
	);
#undef ASSIGN
}

/* CASE statement (load-case) */

cbool parseGRAV(parse *p, case_stmt *c){
	vec3 g;
	return (
		parseComments(p)
		&& parseThisString(p,"GRAV")/* && assign(fprintf(stderr,"GRAV")) */
		&& parseWS(p) && assign(fprintf(stderr," "))
		&& parseDouble(p,&g.x) //&& assign(fprintf(stderr,"%d",id))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,&g.y) //&& assign(fprintf(stderr,"%d",id))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,&g.z) //&& assign(fprintf(stderr,"%d",id))
		&& parseEOLplus(p) //&& assign(fprintf(stderr," "))
		&& case_add_gravity(c,g)
		//&& assign(fprintf(stderr,"%e %e %e\n", gx,gy,gz))
	);
}

cbool parseNDLD(parse *p, case_stmt *c){
	unsigned nodeid;
	vec3 F, M;
	return (
		parseComments(p)
		&& parseThisString(p,"NDLD")// && assign(fprintf(stderr,"NDLD"))
		&& parseWS(p)// && assign(fprintf(stderr," "))
		&& parseNumber(p,&nodeid)
		//&& assign(fprintf(stderr," case=%d node=%d",c->id,nodeid))
		&& parseWS(p)
		//&& assign(fprintf(stderr," Fx=%e Fy=%e Fz=%e Mx=%e My=%e Mz=%e\n",c->id,nodeid,Fx,Fy,Fz,Mx,My,Mz))
		&& parseDouble(p,&F.x)
		&& parseWS(p)
		&& parseDouble(p,&F.y)
		&& parseWS(p)
		&& parseDouble(p,&F.z)
		&& parseWS(p)
		&& parseDouble(p,&M.x)
		&& parseWS(p)
		&& parseDouble(p,&M.y)
		&& parseWS(p)
		&& parseDouble(p,&M.z)
		&& parseEOLplus(p)
		&& case_add_node_load(c,nodeid, F, M)
		//&& assign(fprintf(stderr,"Added node %d to case %d, now has %d loads\n",nodeid,c->id, ARRAY_NUM(c->data)))
	);
}

/*
	Two possible member load statements that we need to parse and ignore:s
		MBLD  6000  CONC  MY  GL  FR         6.0000       0.5000
		MBLD  6000  CONC  FZ  GL  FR        -3.0000       0.5000
*/
cbool parseMBLD(parse *p, case_stmt *c){
	unsigned membid;
	const char FZ[] = "FZ";
	const char MY[] = "MY";
	const char *dir;
	double a, b;
	static cbool warnignored = 0;
#define ASSIGN(Z) 1
	return (
		parseComments(p)
		&& parseThisString(p,"MBLD") && ASSIGN(fprintf(stderr,"MBLD"))	
		&& parseWS(p) && ASSIGN(fprintf(stderr,"NDLD"))
		&& parseNumber(p,&membid)
		&& parseWS(p) 
		&& (
			parseThisString(p,"CONC")
		) && parseWS(p)  && (
			(parseThisString(p,MY) && assign(dir=MY))
			|| (parseThisString(p,FZ) && assign(dir=FZ))
		) && parseWS(p) && (
			parseThisString(p,"GL")
		) && parseWS(p) && (
			parseThisString(p,"FR")
		)
		&& parseWS(p)
		&& parseDouble(p,&a) /* what's this? */
		&& parseWS(p)
		&& parseDouble(p,&b) /* what's this? */
		&& parseEOLplus(p) //&& assign(fprintf(stderr," "))
		&& (warnignored || (
			assign(fprintf(stderr,"WARNING: 'MBLD' statement(s) ignored.\n")) 
			&& assign(warnignored=1)
		))
#undef ASSIGN
	);
}


cbool parseCOMB(parse *p, model *a, case_stmt *c){
	unsigned subcaseid;
	double factor;
#define ASSIGN(Z) 1
	return ((
		parseComments(p) && ASSIGN(fprintf(stderr,"Adding to case %d\n",c->id))
		&& parseThisString(p,"COMB") && ASSIGN(fprintf(stderr,"COMB"))
	) && (
		(
			parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseNumber(p,&subcaseid) && ASSIGN(fprintf(stderr,"%d",subcaseid))
			&& parseWS(p) && ASSIGN(fprintf(stderr," "))
			&& parseDouble(p,&factor) && ASSIGN(fprintf(stderr,"%f",factor))
			&& parseEOLplus(p) && ASSIGN(fprintf(stderr,"\n"))
			&& case_add_comb(a,c,subcaseid,factor)
		) || (
			assign(fprintf(stderr,"Failed to parse COMB\n"))
		)
	));
#undef ASSIGN
}

cbool parseCASE(parse *p, model *a){
	case_stmt c;
	c.type = CASE_UNDEFINED;
	//ndld_stmt ndld[MAXNDLDS];
#define ASSIGN(Z) 1
	return (
		parseComments(p)
		&& parseThisString(p,"CASE") && ASSIGN(fprintf(stderr,"CASE"))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseNumber(p,&c.id) && ASSIGN(fprintf(stderr,"%d",c.id))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseStrExcept(p,"\n\r",c.name,MAXCASENAME) && ASSIGN(fprintf(stderr,"\"%s\"",c.name))
		&& parseEOLplus(p) && ASSIGN(fprintf(stderr,"\n"))
		&& (
			(
				parseGRAV(p,&c)
				&& many((parseNDLD(p,&c) || parseMBLD(p,&c)))
			) || (
				one_or_more((parseNDLD(p,&c) || parseMBLD(p,&c)))
			) || (
				one_or_more(parseCOMB(p,a,&c))
			)
		)
		&& (
			model_add_case(a,&c)
			|| assign(fprintf(stderr,"Unable to add case to model\n"))
		)
	);
#undef ASSIGN
}

/* MICROSTRAN ARCHIVE FILE*/

cbool parseModelMicrostran(parse *p, model **a){
	unsigned version, type, vert;
	unit_stmt unit;
	model *a1 = NULL;

	if(!p)return fail;

#define ASSIGN(Z) 1
	return (
		parseNumericLine(p,"VERS",&version) //&& assign(fprintf(stderr,"VERS %d\n",version))
		&& parseNumericLine(p,"TYPE",&type) //&& assign(fprintf(stderr,"TYPE %d\n",type))
		&& parseNumericLine(p,"VERT",&vert) //&& assign(fprintf(stderr,"VERT %d\n",vert))
		&& parseUNIT(p,&unit) //&& assign(fprintf(stderr,"UNIT parsed\n"))
		&& assign(a1 = model_create(version, type, vert, unit))
		//&& assign(fprintf(stderr,"model created...\n"))
		&& one_or_more(parseNODE(p,a1)) && ASSIGN(fprintf(stderr,"Parsed %d NODEs\n",a1->num_nodes))
		&& one_or_more(parseMEMB(p,a1)) && ASSIGN(fprintf(stderr,"Parsed %d MEMBs\n",a1->num_membs))
		&& many(parseMOFF(p,a1)) && ASSIGN(fprintf(stderr,"MOFFs parsed\n"))
		&& many(parseMSPR(p,a1)) && ASSIGN(fprintf(stderr,"MSPRs parsed\n"))
		&& many(parseMTYP(p,a1)) && ASSIGN(fprintf(stderr,"MTYPs parsed\n"))
		&& one_or_more(parsePROP(p,a1)) && ASSIGN(fprintf(stderr,"Parsed %d PROPs\n",a1->num_props))
		&& many(parseMATL(p,a1)) && ASSIGN(fprintf(stderr,"MATLs parsed\n"))
		&& many(parseCASE(p,a1)) && ASSIGN(fprintf(stderr,"CASEs parsed\n"))
		&& parseThisString(p,"END") && ASSIGN(fprintf(stderr,"END parsed\n"))
		&& parseEOLplus(p) && ASSIGN(fprintf(stderr,"\nEND parsed\n"))
		&& assign(*a = a1)
	) || ((a1?model_destroy(a1):0), fail);
#undef ASSIGN
}


