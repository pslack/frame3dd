#define MSTRANP_BUILD
#include "case.h"
#include "model.h"

#include "new.h"
#include <string.h>
#include <stdio.h>

/* CASE methods */

case_stmt *case_create(unsigned caseid, const char *name){
	case_stmt *c;
	c = NEW(case_stmt);
	c->id = caseid;
	c->type = CASE_UNDEFINED;
	strncpy(c->name, name, MAXCASENAME);
	if(strlen(name)>= MAXCASENAME){
		c->name[MAXCASENAME-1]='\0';
	}
	return c;
}

case_stmt *case_copy(const case_stmt *oc){
	case_stmt *c;
	c = NEW(case_stmt);
	c->id = oc->id;
	c->type = oc->type;
	strncpy(c->name, oc->name, MAXCASENAME);
	if(strlen(oc->name)>= MAXCASENAME){
		c->name[MAXCASENAME-1]='\0';
	}
	if(c->type != CASE_UNDEFINED){
		c->data = array_copy(oc->data);
	}
	return c;
}

cbool case_add_comb(
	model *a, case_stmt *c, unsigned subcaseid, double factor
){
	comb_stmt co;
	if(c->type == CASE_LOADS){
		fprintf(stderr,"case %d is of type LOADS, can't include COMB declaration\n",c->id);
		return 0;
	}
	if(c->type == CASE_UNDEFINED){
		c->type = CASE_COMB;
		c->data = ARRAY_CREATE(comb_stmt,5);
	}

	co.c = model_find_case(a, subcaseid);
	if(!co.c){
		fprintf(stderr,"Subcase %d not found in COMB case %d\n",subcaseid,c->id);
		return 0;
	}

	co.factor = factor;

	if(!array_append(&(c->data), &co)){
		fprintf(stderr,"failed to append subcase to case %d\n",c->id);
		return 0;
	}

	return 1;
}

cbool case_add_gravity(case_stmt *c,double gx, double gy, double gz){
	if(c->type == CASE_COMB){
		fprintf(stderr,"case %d is of type COMB, can't include GRAV declaration\n",c->id);
		return 0;
	}
	if(c->type == CASE_UNDEFINED){
		c->type = CASE_LOADS;
		c->data = ARRAY_CREATE(ndld_stmt,5);
	}

	c->gx = gx;
	c->gy = gy;
	c->gz = gz;

	return 1;
}


cbool case_add_node_load(case_stmt *c,unsigned nodeid
	, double Fx, double Fy, double Fz, double Mx, double My, double Mz){
	ndld_stmt n;

	if(c->type == CASE_COMB){
		fprintf(stderr,"case %d is of type COMB, can't include NDLD declaration\n",c->id);
		return 0;
	}
	if(c->type == CASE_UNDEFINED){
		c->type = CASE_LOADS;
		c->data = ARRAY_CREATE(ndld_stmt,5);
	}

	if(case_find_node(c,nodeid)){
		fprintf(stderr,"node %d already present in case %d\n",nodeid,c->id);
		return 0;
	}

	n.Fx = Fx;
	n.Fy = Fy;
	n.Fz = Fz;
	n.Mx = Mx;
	n.My = My;
	n.Mz = Mz;
	n.node = nodeid;

	//fprintf(stderr,"Node %d being appended...\n",n.node);

	if(!array_append(&(c->data), &n)){
		fprintf(stderr,"Failed to append node %d...\n",n.node);
		return 0;
	}
	return 1;
}

ndld_stmt *case_find_node(case_stmt *c, unsigned nodeid){
	unsigned i;
	if(c->type != CASE_LOADS){
		fprintf(stderr,"CASE is not defined yet\n");
		return NULL;
	}
	for(i=0; i < ARRAY_NUM(c->data); ++i){
		ndld_stmt *n = array_get(&(c->data), i);
		if(n->node == nodeid)return n;
	}
	return NULL;
}


