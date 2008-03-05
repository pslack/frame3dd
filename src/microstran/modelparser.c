/** @FILE
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
	unsigned nodeid, flags;
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
		&& parseNumber(p,&flags)
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
	/* FIXME currently these are not being recorded into the 'model' data structure */
	unsigned memberid;
	vec3 deltaA, deltaB;
	const char *code;
	const char *LO = "LO";
	static int warnMOFF = 0;
	return (
		parseComments(p)
		&& parseThisString(p,"MOFF")
		&& parseWS(p)
		&& parseNumber(p,&memberid)
		&& parseWS(p)
		&& parseThisString(p,"LO") && assign(code=LO)
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
		&& model_add_member_offset(a, memberid, code, deltaA, deltaB)
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
		&& parseThisString(p,"TONLY")
		&& parseEOLplus(p)
		&& assign(warnMTYP ? 0 : (fprintf(stderr,"WARNING: 'MTYP' statements ignored.\n"),warnMTYP=1))
	);
}	

/* PROP statement */

cbool parsePROP(parse *p, model *a){
	unsigned id;
	char libr[MAXPROPLIBNAME] = "";
	char name[MAXPROPNAME] = "";
	char desc[MAXPROPDESC] = "";
	cbool isdefault;
	double vals[6];
	return (
		parseComments(p)
		&& parseThisString(p,"PROP") //&& assign(fprintf(stderr,"PROP"))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&id) //&& assign(fprintf(stderr,"%d",id))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& ((
			parseThisString(p,"LIBR") //&& assign(fprintf(stderr,"LIBR"))
			&& parseWS(p) //&& assign(fprintf(stderr," "))
			&& parseStrExcept(p," \n\r\t",libr,MAXPROPLIBNAME) //&& assign(fprintf(stderr,"%s",libr))
			&& parseWS(p) //&& assign(fprintf(stderr," "))
			&& parseStrExcept(p," \n\r\t",name,MAXPROPNAME) //&& assign(fprintf(stderr,"%s",name))
			&& parseWS(p) //&& assign(fprintf(stderr," "))
			&& parseThisChar(p,'Y') //&& assign(fprintf(stderr,"Y"))
			&& maybe(parseWS(p) && parseThisString(p,"default") && assign(isdefault=1)/*&& assign(fprintf(stderr," default"))*/) 
			&& parseEOLplus(p) //&& assign(fprintf(stderr,"\n\t"))
		) || (
			parseThisString(p,"PRIS") //&& assign(fprintf(stderr,"PRIS"))
			&& parseWS(p) //&& assign(fprintf(stderr," "))
			&& parseStrExcept(p," \n\r\t",name,MAXPROPNAME) //&& assign(fprintf(stderr,"%s",name))
			&& parseWS(p) //&& assign(fprintf(stderr," "))
			&& parseStrExcept(p,"\n\r\t",desc,MAXPROPDESC) //&& assign(fprintf(stderr,"%s",desc))
			&& parseEOLplus(p)
		)) //&& assign(fprintf(stderr,"\n\t"))
		&& parseDouble(p,vals+0) //&& assign(fprintf(stderr,"%e",vals[0]))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,vals+1) //&& assign(fprintf(stderr,"%e",vals[1]))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,vals+2) //&& assign(fprintf(stderr,"%e",vals[2]))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,vals+3) //&& assign(fprintf(stderr,"%e",vals[3]))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,vals+4) //&& assign(fprintf(stderr,"%e",vals[4]))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,vals+5) //&& assign(fprintf(stderr,"%e",vals[5]))
		&& parseEOLplus(p) //&& assign(fprintf(stderr,"\n"))
		&& model_add_prop(a, id, libr, name, desc, isdefault, vals)
	);	
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
	double gx,gy,gz;
	return (
		parseComments(p)
		&& parseThisString(p,"GRAV")/* && assign(fprintf(stderr,"GRAV")) */
		&& parseWS(p) && assign(fprintf(stderr," "))
		&& parseDouble(p,&gx) //&& assign(fprintf(stderr,"%d",id))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,&gy) //&& assign(fprintf(stderr,"%d",id))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,&gz) //&& assign(fprintf(stderr,"%d",id))
		&& parseEOLplus(p) //&& assign(fprintf(stderr," "))
		&& case_add_gravity(c,gx, gy, gz)
		//&& assign(fprintf(stderr,"%e %e %e\n", gx,gy,gz))
	);
}

cbool parseNDLD(parse *p, case_stmt *c){
	unsigned nodeid;
	double Fx,Fy,Fz, Mx,My,Mz;
	return (
		parseComments(p)
		&& parseThisString(p,"NDLD")// && assign(fprintf(stderr,"NDLD"))
		&& parseWS(p)// && assign(fprintf(stderr," "))
		&& parseNumber(p,&nodeid)
		//&& assign(fprintf(stderr," case=%d node=%d",c->id,nodeid))
		&& parseWS(p)
		//&& assign(fprintf(stderr," Fx=%e Fy=%e Fz=%e Mx=%e My=%e Mz=%e\n",c->id,nodeid,Fx,Fy,Fz,Mx,My,Mz))
		&& parseDouble(p,&Fx)
		&& parseWS(p)
		&& parseDouble(p,&Fy)
		&& parseWS(p)
		&& parseDouble(p,&Fz)
		&& parseWS(p)
		&& parseDouble(p,&Mx)
		&& parseWS(p)
		&& parseDouble(p,&My)
		&& parseWS(p)
		&& parseDouble(p,&Mz)
		&& parseEOLplus(p)
		&& case_add_node_load(c,nodeid,Fx,Fy,Fz,Mx,My,Mz)
		//&& assign(fprintf(stderr,"Added node %d to case %d, now has %d loads\n",nodeid,c->id, ARRAY_NUM(c->data)))
	);
}

cbool parseCOMB(parse *p, model *a, case_stmt *c){
	unsigned subcaseid;
	double factor;
	return (
		parseComments(p) //&& assign(fprintf(stderr,"Adding to case %d\n",c->id))
		&& parseThisString(p,"COMB") //&& assign(fprintf(stderr,"COMB"))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseNumber(p,&subcaseid)// && assign(fprintf(stderr,"%d",subcaseid))
		&& parseWS(p) //&& assign(fprintf(stderr," "))
		&& parseDouble(p,&factor) //&& assign(fprintf(stderr,"%f",factor))
		&& parseEOLplus(p) //&& assign(fprintf(stderr,"\n"))
		&& case_add_comb(a,c,subcaseid,factor)
	);
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
		&& parseNumber(p,&(c.id)) && ASSIGN(fprintf(stderr,"%d",c.id))
		&& parseWS(p) && ASSIGN(fprintf(stderr," "))
		&& parseStrExcept(p,"\n\r",c.name,MAXCASENAME) && ASSIGN(fprintf(stderr,"\"%s\"",c.name))
		&& parseEOLplus(p) && ASSIGN(fprintf(stderr,"\n"))
		&& ((maybe(parseGRAV(p,&c))
			&& parseNDLD(p,&c) && many(parseNDLD(p,&c))
		) || (
			parseCOMB(p,a,&c) && many(parseCOMB(p,a,&c))
		))
		&& model_add_case(a,&c)
	);
#undef ASSIGN
}

/* MICROSTRAN ARCHIVE FILE*/

cbool parseModelMicrostran(parse *p, model **a){
	unsigned version, type, vert;
	unit_stmt unit;
	model *a1 = NULL;
#define ASSIGN(Z) 1
	return (
		parseNumericLine(p,"VERS",&version) //&& assign(fprintf(stderr,"VERS %d\n",version))
		&& parseNumericLine(p,"TYPE",&type) //&& assign(fprintf(stderr,"TYPE %d\n",type))
		&& parseNumericLine(p,"VERT",&vert) //&& assign(fprintf(stderr,"VERT %d\n",vert))
		&& parseUNIT(p,&unit) //&& assign(fprintf(stderr,"UNIT parsed\n"))
		&& assign(a1 = model_create(version, type, vert, unit))
		//&& assign(fprintf(stderr,"model created...\n"))
		&& many(parseNODE(p,a1)) && ASSIGN(fprintf(stderr,"NODEs parsed\n"))
		&& many(parseMEMB(p,a1)) && ASSIGN(fprintf(stderr,"MEMBs parsed\n"))
		&& many(parseMOFF(p,a1)) && ASSIGN(fprintf(stderr,"MOFFs parsed\n"))
		&& many(parseMSPR(p,a1)) && ASSIGN(fprintf(stderr,"MSPRs parsed\n"))
		&& many(parseMTYP(p,a1)) && ASSIGN(fprintf(stderr,"MTYPs parser\n"))
		&& many(parsePROP(p,a1)) && ASSIGN(fprintf(stderr,"PROPs parsed\n"))
		&& many(parseMATL(p,a1)) && ASSIGN(fprintf(stderr,"MATLs parsed\n"))
		&& many(parseCASE(p,a1)) && ASSIGN(fprintf(stderr,"CASEs parsed\n"))
		&& parseThisString(p,"END") && parseEOLplus(p) && ASSIGN(fprintf(stderr,"\nEND parsed\n"))
		&& assign(*a = a1)
	) || ((a1?model_destroy(a1):0), fail);
#undef ASSIGN
}


