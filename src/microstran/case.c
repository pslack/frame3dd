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
	c->data = array_create(0,0);
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
	co.factor = factor;

	if(c->type == CASE_LOADS){
		fprintf(stderr,"case %d is of type LOADS, can't include COMB declaration\n",c->id);
		return 0;
	}
	if(c->type == CASE_UNDEFINED){
		c->type = CASE_COMB;
		c->data = ARRAY_CREATE(comb_stmt,5);
	}

	//fprintf(stderr,"CASE %d contains %d COMBs\n",c->id, ARRAY_NUM(c->data));

	if(!model_find_case(a, subcaseid, &co.caseindex)){
		fprintf(stderr,"Subcase %d not found in COMB case %d\n",subcaseid,c->id);
		return 0;
	}

	if(!array_append(&c->data, &co)){
		fprintf(stderr,"failed to append subcase to case %d\n",c->id);
		return 0;
	}

	comb_stmt *co1 = array_get(&c->data, ARRAY_NUM(c->data)-1);
	case_stmt *case1;
	case1 = array_get(&a->cases, co1->caseindex);
	//fprintf(stderr,"Added subcase %d '%s' to case %d\n",case1->id, case1->name, c->id);

	return 1;
}

cbool case_add_gravity(case_stmt *c, vec3 g){
	if(c->type == CASE_COMB){
		fprintf(stderr,"case %d is of type COMB, can't include GRAV declaration\n",c->id);
		return 0;
	}
	if(c->type == CASE_UNDEFINED){
		c->type = CASE_LOADS;
		c->data = ARRAY_CREATE(ndld_stmt,5);
	}

	c->g = g;

	return 1;
}


cbool case_add_node_load(case_stmt *c,unsigned nodeid, vec3 F, vec3 M){
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

	n.node = nodeid;
	n.F = F;
	n.M = M;

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

/*
	in this code, we will assume that each case may contain multiple
	loads on the same node (although in our CASE parser we limited it to one)
*/
cbool case_total_load_node(model *a, case_stmt *c, unsigned nodeid, ndld_stmt *nl){
	ndld_stmt n;
	n.node = nodeid;
	n.F = VEC3_ZERO;
	n.M = VEC3_ZERO;
	unsigned i;
	cbool found = 0;
	/* fprintf(stderr,"Looking for loads on node %d in case %d...\n",nodeid,c->id); */
	switch(c->type){
		case CASE_LOADS:
			/* fprintf(stderr,"Case %d has %d node loads\n",c->id,ARRAY_NUM(c->data)); */
			for(i=0; i<ARRAY_NUM(c->data); ++i){
				ndld_stmt *nl1 = array_get(&c->data, i);
				if(nl1->node == nodeid){
					//VEC3_PR(nl1->F);
					++found;
					n.F = vec3_add(n.F, nl1->F);
					n.M = vec3_add(n.M, nl1->M);
				}
			}
			break;
		case CASE_COMB:
			/* fprintf(stderr,"COMB case with %d subcases\n",ARRAY_NUM(c->data)); */
			for(i=0; i<ARRAY_NUM(c->data); ++i){
				comb_stmt *comb1 = (comb_stmt *)array_get(&(c->data), i);
				case_stmt *case2 = array_get(&a->cases, comb1->caseindex);
				/* fprintf(stderr,"Descending into sub-case %d '%s' (factor = %f)...\n",case2->id,case2->name,comb1->factor); */
				ndld_stmt nl1;
				if(case_total_load_node(a, case2, nodeid, &nl1)){
					++found;
					n.F = vec3_add(n.F, vec3_scale(nl1.F, comb1->factor));
					n.M = vec3_add(n.M, vec3_scale(nl1.M, comb1->factor));
				}
			}
			break;
		case CASE_UNDEFINED:
			break;
	}
	*nl = n;
	//if(found)fprintf(stderr,"%d NDLDs found for node %d\n",found,nodeid);
	return found;
}
		
		

