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

typedef struct node_stmt_{
	unsigned id;
	vec3 pos;
	unsigned flags;
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



typedef struct memb_stmt_{
	unsigned id;
	unsigned fromnode; /**< NOTE: we store index to the 'node' array, rather than the node ID here */
	unsigned tonode; /**< NOTE: we store index to the 'node' array, rather than the node ID here */
	member_orientation orient;
	unsigned prop; /**< NOTE: this is a property ID (direct from microstran) */
	unsigned matl; /**< NOTE: this is a material ID (direct from microstran) */
	unsigned flags1;
	unsigned flags2;
} memb_stmt;

/* MOFF statement (member offsets) */

typedef struct moff_stmt_{
	unsigned id;
	const char *code; /**< not sure what this one means at this stage, only value seen has been 'LO' */
	vec3 deltafrom; /**< offset at the 'from' node end */
	vec3 deltato; /**< offset at the 'to' node end */
} moff_stmt;

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

enum prop_vals{
	PROP_A, PROP_2, PROP_3, PROP_J, PROP_IYY, PROP_IXX,MAXPROPVALS
};

/* PROP statement (section properties) */

typedef struct prop_stmt_{
	unsigned id;
	char libr[MAXPROPLIBNAME];
	char name[MAXPROPNAME];
	char desc[MAXPROPDESC];
	cbool isdefault;
	double vals[MAXPROPVALS];
} prop_stmt;

/* MATL statement (material stiffness/yield etc) */

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

/*
	FIXME these structures are fixed-size, which is a bit silly
*/
typedef struct model_{
	int version;
	int type;
	int vert;
	unit_stmt unit;
	node_stmt node[MAXNODES];
	unsigned num_nodes;
	memb_stmt memb[MAXMEMBS];
	unsigned num_membs;
	array moffs;
	prop_stmt prop[MAXPROPS];
	unsigned num_props;
	matl_stmt matl[MAXMATLS];
	unsigned num_matls;
	array cases;
} model;

model *model_create(unsigned version, unsigned type, unsigned vert, unit_stmt unit);
MSTRANP_API void model_destroy(model *a);

cbool model_add_node(model *a, node_stmt node);

MSTRANP_API cbool model_find_node(const model *a, unsigned id, unsigned *index);

MSTRANP_API cbool model_find_memb(const model *a, const unsigned membid, unsigned *membindex);
MSTRANP_API cbool model_find_memb_from_to(const model *a, const unsigned nodeid1, const unsigned nodeid2, unsigned *membindex);

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
cbool model_add_member_offset(model *a, unsigned memberid, const char *code, vec3 deltafrom, vec3 deltato);

/**
	Find whether a member has an offset applied to its end locations.
	@return pointer to the offset struct, or NULL if no offset exists.
*/
MSTRANP_API moff_stmt *model_find_member_offset(const model *a, const unsigned membid);

cbool model_add_prop(model *a, unsigned id, char libr[], char name[], char desc[]
		, cbool isdefault, double vals[MAXPROPVALS]
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

case_stmt *model_find_case(model *a,unsigned caseid);

cbool model_add_case(model *a, case_stmt *c);

struct casedisplacements_;

MSTRANP_API model *model_copy(const model *m);

/* modify a model by applying displacements to each node position */
MSTRANP_API cbool model_apply_displacements(model *a, struct casedisplacements_ *cd);

#ifdef __cplusplus
}
#endif

#endif /* MSTRANP_MODEL_H */

