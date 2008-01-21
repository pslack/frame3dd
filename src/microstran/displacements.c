#define MSTRANP_BUILD

#include "displacements.h"
#include "new.h"
#include "array.h"
#include "error.h"

#include <stdio.h>
#include <string.h>

modeldisplacements *modeldisplacements_create(){
	modeldisplacements *md;
	md = NEW(modeldisplacements);
	strcpy(md->filename,"");
	md->cases = ARRAY_CREATE(casedisplacements,10);
	return md;
}

void modeldisplacements_destroy(modeldisplacements *d){
	/* FIXME destroy contained data too */
	unsigned i;
	for(i=0;i<ARRAY_NUM(d->cases);++i){
		casedisplacements_destroy(array_get(&(d->cases),i));
	}
	array_destroy(&(d->cases));
	free(d);
	d = NULL;
}

cbool modeldisplacements_add_case(modeldisplacements *md, casedisplacements *cd){
	//fprintf(stderr,"Model displacements: contains %d cases\n",ARRAY_NUM(md->cases));
	array_append(&(md->cases),cd);
}

casedisplacements *modeldisplacements_get_case(modeldisplacements *md, unsigned index){
	return array_get(&(md->cases), index);
}

void casedisplacements_destroy(casedisplacements *cd){
	array_destroy(&(cd->nodes));
	free(cd);
	cd = NULL;
}

casedisplacements *casedisplacements_create(unsigned id, const char *name){
	casedisplacements *cd;
	cd = NEW(casedisplacements);
	if(!cd)abort();
	cd->nodes = ARRAY_CREATE(nodedisplacement,10);
	cd->id = id;
	strncpy(cd->name,name,MAXCASENAME);
	cd->name[MAXCASENAME-1]='\0';
	return cd;
}

cbool casedisplacements_add_node(
		casedisplacements *cd, unsigned nodeid
			, double dx, double dy, double dz, double rx, double ry, double rz
){
	nodedisplacement nd;
	//fprintf(stderr,"ADDING NODE %d to case %d",nodeid,cd->id);
	nd.node = nodeid;
	nd.dx = dx;
	nd.dy = dy;
	nd.dz = dz;
	nd.rx = rx;
	nd.ry = ry;
	nd.rz = rz;
	array_append(&(cd->nodes),&nd);
	return 1;
}

casedisplacements *modeldisplacements_find_case(modeldisplacements *md, unsigned caseid){
	unsigned i;
	//fprintf(stderr,"Model displacement contains %d cases\n",ARRAY_NUM(md->cases));
	casedisplacements *cd;
	for(i=0; i<ARRAY_NUM(md->cases); ++i){
		cd = array_get(&(md->cases), i);
		if(cd->id == caseid){
			return cd;
		}
	}
	return NULL;
}

nodedisplacement *casedisplacements_find_node(casedisplacements *cd, const unsigned nodeid){
	unsigned i;
	nodedisplacement *nd;
	for(i=0; i<ARRAY_NUM(cd->nodes); ++i){
		nd = array_get(&(cd->nodes), i);
		if(nd->node == nodeid){
			return nd;
		}
	}
	return NULL;
}


