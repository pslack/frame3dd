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
*//** @file
	Data structures to hold model case deflection data
*/
#ifndef DISPLACEMENTS_H
#define DISPLACEMENTS_H

#include "config.h"
#include "model.h"

#ifdef __cplusplus
extern "C"{
#endif

#define MAXDISPLACEMENTSNAME 80

/**
	Displacement (both translation and rotation)
	of a node, together with the node's ID.
*/
typedef struct nodedisplacement_{
	int node;
	double dx,dy,dz; /* displacements */
	double rx,ry,rz; /* rotations */
} nodedisplacement;

/**
	Displacement data for a specific load case.
*/
typedef struct casedisplacements_{
	int id;
	char name[MAXCASENAME];
	array nodes; /* array (list) of node_displacement */
} casedisplacements;

/**
	Displacement data structure for containing results parsed from
	.p1 data files.
*/
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

