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
*//**
	@file
	specification of model data structures and methods.
*//**
	@page mstranp Microstran parser
	Code in this directory provides the ability to load Microstran .arc model
	files as well as .p1 files containing member forces and node displacements
	as computed by Microstran. Ultimately we aim to use this code to
	interface with the main FRAME3DD calculation engine, so that FRAME3DD can 
	also be used to analyse models generated in the .arc file format. Commercial
	software that supports the .arc format includes Microstran, Spacegass, 
	and Multiframe.

	Key data structures in this section are:
	@li model
    @li modelforces
    @li modeldisplacements
	@li parse
*/

#ifndef MSTRANP_MODEL_H
#define MSTRANP_MODEL_H

#include "config.h"
#include "types.h"
#include "case.h"
#include "array.h"
#include "vec3.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

#define MAXUNIT 10

/* NODE structure */

#define MSTRANP_NODE_FIXX  0x01
#define MSTRANP_NODE_FIXY  0x02
#define MSTRANP_NODE_FIXZ  0x04
#define MSTRANP_NODE_FIXMX 0x08
#define MSTRANP_NODE_FIXMY 0x10
#define MSTRANP_NODE_FIXMZ 0x20

/**
	Node statement structure.
*/
typedef struct node_stmt_{
	unsigned id; /**< node ID (from Microstran .arc) */
	vec3 pos; /**< node position */
	unsigned flags; /**< composed of ORed values like MSTRANP_NODE_FIX*. */
} node_stmt;

MSTRANP_API node_stmt node_create(unsigned nodeid,vec3 pos,unsigned flags);
MSTRANP_API int node_print(FILE *f, const node_stmt *n);

/* changes the node 'n' and also returns the pointer back again */
node_stmt *node_translate(node_stmt *n, double dx, double dy, double dz);

/* MEMB structure */

/**
	Structure to hold data on member orientation. In Microstran, you can specify
	alignment to one of the global coordinate system axes, or else alignment
	to a node.
*/
typedef struct member_orientation_struct{
	char axis; /** Axis 'X', 'Y', or 'Z', or '\0' to indicate orientation to a node */
	union{
		unsigned node; /**< if axis=='\0', this will hode an index into the 'node' array (NOT a node ID) */
		char dir; /**< if orientation to an axis, this will hold '+' or '-' for the direction of orientation */
	};
} member_orientation;


/**
	Data structure for frame members imported from Microstran.
*/
typedef struct memb_stmt_{
	unsigned id; /**< Member ID */
	unsigned fromnode; /**< NOTE: we store index to the 'node' array, rather than the node ID here */
	unsigned tonode; /**< NOTE: we store index to the 'node' array, rather than the node ID here */
	member_orientation orient; /**< Member orientation */
	unsigned prop; /**< NOTE: this is a property ID (direct from microstran) */
	unsigned matl; /**< NOTE: this is a material ID (direct from microstran) */
	unsigned flags1; /**< flags for the 'from' end? */
	unsigned flags2; /**< flags for the 'to' end? */
} memb_stmt;

/*----------------------
	MOFF statement (member offsets) 
*/

/**
	Enumeration of the possible coordinate systems used in member offsets
*/
typedef enum{
	MSTRANP_COORDS_LOCAL=0 /**< local coordinate system */
	,MSTRANP_COORDS_GLOBAL /**< global coordinate system */
} coord_sys_t;

/**
	Member offset statement, for members whose ends do not finish exactly
	at a node location. These structures are held separately from the member
	data, because not all members have an offset, and because of the way that
	data is read in from the Microstran .arc file.
*/
typedef struct moff_stmt_{
	unsigned id;
	coord_sys_t coordsys; /**< coordinate system for these offsets, either local or global */
	vec3 deltafrom; /**< offset at the 'from' node end */
	vec3 deltato; /**< offset at the 'to' node end */
} moff_stmt;

MSTRANP_API int moff_print(FILE *f, const moff_stmt *o);

/**
	UNIT statement, which specifies the units of measurement in the
	Microstran .arc file. At present we don't make use of these data in
	subsequent output.
*/
typedef struct unit_stmt_{
	unsigned num;
	char lengthunit[MAXUNIT];
	char forceunit[MAXUNIT];
	char massunit[MAXUNIT];
	char tempunit[MAXUNIT];
} unit_stmt;

#define MAXPROPLIBNAME 10
#define MAXPROPNAME 40
#define MAXPROPDESC 200

/*
	Enumeration of the steel properties that are kept in the prop_stmt
	'vals' array.
*/
enum prop_vals{
	PROP_A, PROP_2, PROP_3, PROP_J, PROP_IYY, PROP_IXX,MAXPROPVALS
};

/**
	PROP statement (section properties)
*/
typedef struct prop_stmt_{
	unsigned id;
	char libr[MAXPROPLIBNAME];
	char name[MAXPROPNAME];
	char desc[MAXPROPDESC];
	double vals[MAXPROPVALS];
} prop_stmt;

/**
	MATL statement (material stiffness/yield etc)
*/
typedef struct matl_stmt_{
	unsigned id;
	double E; /* force unit / length_unit^2 */
	double sigma_y; /* force unit / length_unit^2 */
	double rho; /* mass_unit / length_unit^3 */
	double beta; /* temp_unit^-1 ...guessing about this one */
} matl_stmt;

/* Overall Microstran 'model' data type */

#define MAXNODES 10000
#define MAXMEMBS 10000
#define MAXPROPS 500
#define MAXMATLS 50
#define MAXCASES 200

/**
	Top-level data structure for Microstran model data.

	@TODO FIXME these structures are fixed-size, which is a bit silly
*/
typedef struct model_{
	int version; /**< file format, currently expected to be '5' */
	int type; /**< type of model, currently expected to be '5' */
	int vert; /**< number of spatial dimensions? Expected to be '3' */
	unit_stmt unit; /**< units of measurement for this model */
	node_stmt node[MAXNODES];
	unsigned num_nodes; /**< number of nodes in the model */
	memb_stmt memb[MAXMEMBS];
	unsigned num_membs; /**< number of members in the model */
	array moffs; /**< dynamically-sized array of member offsets (of type moff_stmt) */
	prop_stmt prop[MAXPROPS]; /** section properties (steel cross section data) */
	unsigned num_props; /**< number of cross-sections used in the model */
	matl_stmt matl[MAXMATLS];
	unsigned num_matls; /**< number of materials (steel properties etc */
	array cases; /**< dynamically-sized array of load case data (of type case_stmt) */
} model;

model *model_create(unsigned version, unsigned type, unsigned vert, unit_stmt unit);
MSTRANP_API void model_destroy(model *a);

cbool model_add_node(model *a, node_stmt node);

MSTRANP_API cbool model_find_node(const model *a, unsigned id, unsigned *index);

MSTRANP_API cbool model_find_memb(const model *a, const unsigned membid, unsigned *membindex);

/**
	Find a member that bridges between nodeid1 and nodeid2. Flips the 'from' and 'to' if necessary, in order to find a match.
	The behaviour for the case where a model contains multiple members between a given node is undefined.

	@param nodeid1 'from' node ID
	@param nodeid2 'to' node ID
	@param membindex (returned) index into the model->memb array (NOTE: NOT equal to the member ID)
	@return 1 if such a member is found (in which case, *membindex will have been set), 0 otherwise.
*/
MSTRANP_API cbool model_find_memb_from_to(const model *a, const unsigned nodeid1, const unsigned nodeid2, unsigned *membindex);

/**
	Find the next member that has either 'from' or 'to' node passing through
	node 'nodeid1', starting at the member pointed to by 'start' (which must
	be a pointer inside the array a->memb. If 'start' is NULL, assume that we
	must start 'just before' a->memb[0].

	This function can be used to (fairly) efficiently iterate through all of
	the members	connected to a node.

	@return pointer into a->memb of next member AFTER start that is connected
	to node 'nodeid', or NULL if no such node exists.
	@param nodeid the node ID which we want the connected members for.
	@param start pointer to the last found member, or NULL if we should start
		at the start.
*/
MSTRANP_API memb_stmt *model_find_memb_from(model *a, const unsigned nodeid1, memb_stmt *start);

cbool model_add_memb(model *a, unsigned id,unsigned fromnode
		,unsigned tonode, member_orientation orient, unsigned prop, unsigned matl
		,unsigned flags1,unsigned flags2
);

/**
	Add member offsets in the form of variations in the positions of the
	member ends relative to the nodes at which the member ends are anchored.
	We're not sure yet how Microstran calculates with these offsets, but
	we want to be able to parse them so that we can correctly _render_ frames that
	include such offsets.
*/
cbool model_add_member_offset(model *a, unsigned memberid, coord_sys_t code, vec3 deltafrom, vec3 deltato);

/**
	Find whether a member has an offset applied to its end locations.
	@return pointer to the offset struct, or NULL if no offset exists.
*/
MSTRANP_API moff_stmt *model_find_member_offset(const model *a, const unsigned membid);

#if 0
/**
	Calculate member offsets in global coordinates (returned in *moff). Useful
	in rendering code.
	@return 1 on success.
*/
MSTRANP_API cbool model_get_member_offset_global(const model *a, const unsigned memberid, moff_stmt *moff);
#endif

cbool model_add_prop(model *a, unsigned id, char libr[], char name[], char desc[]
		, double vals[MAXPROPVALS]
);

/**
	Return the direction of the 'X' (minor?) axis of the member as specified
	by its orientation data. The returned direction vector will be normal
	to the axis of the member. Note that section outline is described in these
	'X' and 'Y' coordinates (see sections.h)
*/
MSTRANP_API vec3 memb_get_orientation(const model *a, const memb_stmt *m);

MSTRANP_API prop_stmt *model_find_prop(model *a, unsigned propid);

cbool model_add_matl(model *a, unsigned id, double E, double sigma_y
		, double rho, double beta
);

MSTRANP_API void model_write_inventory(model *a);

MSTRANP_API cbool model_find_case(model *a,unsigned caseid, unsigned *caseindex);

/**
	Get case at position 'caseindex'. Note that this is NOT the case ID, it
	is the position in the case array instead, for ease of iterating through
	the available cases.

	@note that case pointers are not persistent; if more cases are added to the 
	model, it is possible that the case_stmt* may be invalidated. @endnote

	@return pointer to the case data structure, or NULL if the caseindex is invalid.
*/
MSTRANP_API case_stmt *model_get_case(model *a, unsigned caseindex);

cbool model_add_case(model *a, case_stmt *c);

struct casedisplacements_;

MSTRANP_API model *model_copy(const model *m);

/* modify a model by applying displacements to each node position */
MSTRANP_API cbool model_apply_displacements(model *a, struct casedisplacements_ *cd);

#ifdef __cplusplus
}
#endif

#endif /* MSTRANP_MODEL_H */

