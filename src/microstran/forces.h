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
	Data structures to hold member-end force data
*/
#ifndef DISPLACEMENTS_H
#define DISPLACEMENTS_H

#include "config.h"
#include "model.h"
#include "vec3.h"

#ifdef __cplusplus
extern "C"{
#endif

#define MAXFORCESNAME 80

/**
	Forces and moments on a certain member due to the specified node.
	Each member will have two such force struct, one for each end. Forces
	are in local member coordinate systems, as can be accessed using the
	function 'memb_get_orientation'.

	@NOTE carefully that when a membernodeforce corresponds to a 'from' node,
	the true direction of the force and moment is actually the negative
	of that stored here, because of the sign convention used by Microstran
	for member forces.

	This comes about because positive axial forces are those which cause a
	member to be in tension, regardless of whether one is at the +Z or -Z 
	end of the member. @ENDNOTE

	@NOTE ALSO that in libmicrostranparser's sign convention, local coordinates 
	are 'Z' aligned with the member axis in the direction from-->to, and
	'Y' is the major axis of the cross-section. 'Z' and 'X' are reversed
	compared with the Microstran sign convention. @ENDNOTE
*/
typedef struct membernodeforce_{
	unsigned node;
	vec3 F; /**< force applied on member due to this node */
	vec3 M; /**< moment applied on member due to this node */
} membernodeforce;

/**
	Structure containing all forces calculated for a particular member.
	There are forces and moments at each end of the member.
*/
typedef struct memberforce_{
	unsigned member; /* member ID, as read from the .p1 file */
	membernodeforce fromnode; /* forces for the first-listed node (FIXME check that this node corresponds to the 'from' node in the model struct) */
	membernodeforce tonode; /* forces for the second-listed node (FIXME check that this node corresponds to the 'to' node in the model struct) */
} memberforce;

/**
	All forces from a certain load case. Note that Microstran sometimes doesn't
	export forces for ALL members in a model, so we can't assume that all
	forces will be present for all nodes/members in subsequent processing.
*/
typedef struct caseforces_{
	unsigned id; /**< case ID, as read from the .p1 file. */
	char name[MAXCASENAME]; /**< load case name, as read from the .p1 file. */
	array members; /**< array (list) of memberforce */
} caseforces;

/**
	All forces for all load cases for a particular model. We'll expect to be
	reading this data from a Microstran output '.p1' file, so we also store
	a filename here in this structure.
*/
typedef struct modelforces_{
	char filename[MAXFORCESNAME]; /**< filename from which the data was read. */
	array cases; /* array (list) of caseforces */
} modelforces;

caseforces *caseforces_create(unsigned id, const char *name);
void caseforces_destroy(caseforces *cd);

modelforces *modelforces_create();
void modelforces_destroy(modelforces *mf);

membernodeforce membernodeforce_create(unsigned nodeid, vec3 F, vec3 M);

cbool caseforces_add_member(caseforces *cf, unsigned memberid, membernodeforce fromnode, membernodeforce tonode);

cbool modelforces_add_case(modelforces *mf, caseforces *cf);

/* QUERY FUNCTIONS */

MSTRANP_API caseforces *modelforces_find_case(modelforces *mf, unsigned caseid);

MSTRANP_API memberforce *caseforces_find_member(caseforces *cf, const unsigned memberid);

MSTRANP_API caseforces *modelforces_get_case(modelforces *mf, unsigned index);

#ifdef __cplusplus
};
#endif

#endif /* FORCES_H */
