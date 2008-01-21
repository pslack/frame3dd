/** @FILE
	data structures to hold model case deflection data
*/
#ifndef DISPLACEMENTS_H
#define DISPLACEMENTS_H

#include "config.h"
#include "model.h"

#ifdef __cplusplus
extern "C"{
#endif

#define MAXDISPLACEMENTSNAME 80

typedef struct nodedisplacement_{
	int node;
	double dx,dy,dz; /* displacements */
	double rx,ry,rz; /* rotations */
} nodedisplacement;

typedef struct casedisplacements_{
	int id;
	char name[MAXCASENAME];
	array nodes; /* array (list) of node_displacement */
} casedisplacements;

typedef struct modeldisplacements_{
	char filename[MAXDISPLACEMENTSNAME];
	array cases; /* array (list) of displacement_case */
} modeldisplacements;

casedisplacements *casedisplacements_create(unsigned id, const char *name);
void casedisplacements_destroy(casedisplacements *cd);

cbool casedisplacements_add_node(casedisplacements *cd, unsigned nodeid, double dx, double dy, double dz, double rx, double ry, double rz);

modeldisplacements *modeldisplacements_create();
void modeldisplacements_destroy(modeldisplacements *a);

cbool modeldisplacements_add_case(modeldisplacements *md, casedisplacements *cd);

MSTRANP_API casedisplacements *modeldisplacements_find_case(modeldisplacements *md, unsigned caseid);

nodedisplacement *casedisplacements_find_node(casedisplacements *cd, const unsigned nodeid);

MSTRANP_API casedisplacements *modeldisplacements_get_case(modeldisplacements *md, unsigned index);

#ifdef __cplusplus
}
#endif

#endif

