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
	Load cases data structures, and associated methods.
*/
#ifndef MSTRANP_CASE_H
#define MSTRANP_CASE_H

#include "config.h"
#include "types.h"
#include "array.h"
#include "vec3.h"

#ifdef __cplusplus
extern "C"{
#endif

/* Load case statements */

#define MAXNDLDS 10000

/**
	Node load (NDLD) statement.
*/
typedef struct ndld_stmt_{
	unsigned node;
	vec3 F;
	vec3 M;
} ndld_stmt;

struct case_stmt_;

/**
	Load sub-case specified (for use in COMB statements)
*/
typedef struct comb_stmt_{
	unsigned caseindex; /**< this is the index for the referenced load case
		in the model->cases array. We can't use a pointer directly because
		the array grows dynamically and is relocated when this happens, thus
		invalidating pointers. */
	double factor;
} comb_stmt;

#define MAXCASENAME 200

typedef enum case_type_{
	CASE_UNDEFINED, CASE_LOADS, CASE_COMB
} case_type;

/**
	Structure to store load case, whether it be a 'combination' load case,
	or a 'basic' load case comprised of node loads.

	@TODO FIXME we currently only support NDLD node loads, no support yet
	exists for member loads (MBLD statements in the .arc file).
*/
typedef struct case_stmt_{
	unsigned id;
	char name[MAXCASENAME];
	case_type type;
	unsigned num_sub;
	vec3 g; /**< gravity magnitude/direction, if applicable from GRAV stmt */
	array data;/**< array of ndld_stmt or comb_stmt (depending on case type) */
} case_stmt;

/**
	Find a node-load in a load case, by node ID
*/
ndld_stmt *case_find_node(case_stmt *c, unsigned nodeid);

/**
	Create a new load case, with the specified name
*/
case_stmt *case_create(unsigned caseid, const char *name);

/**
	Duplicate load case data
*/
case_stmt *case_copy(const case_stmt *c);

/**
	Specify the body force (normally gravity) on the structure's mass
*/
cbool case_add_gravity(case_stmt *c, vec3 g);

/**
	Add a node load to a load case
	@param nodeid loaded node ID
	@param F applied force on the node
	@param M applied moment on the node
*/
cbool case_add_node_load(case_stmt *c, unsigned nodeid, vec3 F, vec3 M);

struct model_;

/**
	Add a sub-case to a 'combined' load case. This is used to implement the 
	'COMB' statement from the .arc file format. This allows superposition
	of different load cases for more complex loading scenarios.
*/
cbool case_add_comb(struct model_ *a, case_stmt *c, unsigned subcaseid, double factor);

/**
	Get the total load on a node for a particular node. If there is
	no applied load on the node, return false. Otherwise, insert the total
	applied node load (force, moment) into the ndld_stmt data structure
	pointed to by 'nl'.
	
	@param nodeid node ID of the node we're looking at.
	@param nl (returned) applied force and moment on the node, set to zero
		if no applied load is found.
	@return 1 if node is loaded, 0 otherwise.
*/
MSTRANP_API cbool case_total_load_node(struct model_ *a, case_stmt *c, unsigned nodeid, ndld_stmt *nl);

#ifdef __cplusplus
};
#endif

#endif /* MSTRANP_CASE_H */

