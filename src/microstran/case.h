#ifndef CASE_H
#define CASE_H

#include "config.h"
#include "types.h"
#include "array.h"

/* Load case statements */

#define MAXNDLDS 10000
typedef struct ndld_stmt_{
	unsigned node;
	double Fx,Fy,Fz,Mx,My,Mz;
} ndld_stmt;

struct case_stmt_;

typedef struct comb_stmt_{
	struct case_stmt_ *c;
	double factor;
} comb_stmt;

#define MAXCASENAME 200

typedef enum case_type_{
	CASE_UNDEFINED, CASE_LOADS, CASE_COMB
} case_type;

typedef struct case_stmt_{
	unsigned id;
	char name[MAXCASENAME];
	case_type type;
	unsigned num_sub;
	double gx,gy,gz;
	array data;/* array of ndld_stmt or comb_stmt (depending on case type) */
} case_stmt;

ndld_stmt *case_find_node(case_stmt *c, unsigned nodeid);

case_stmt *case_create(unsigned caseid, const char *name);

case_stmt *case_copy(const case_stmt *c);

cbool case_add_gravity(case_stmt *c,double gx, double gy, double gz);

cbool case_add_node_load(case_stmt *c, unsigned nodeid
	, double Fx, double Fy, double
	, double Mx, double My, double Mz
);

struct model_;

cbool case_add_comb(struct model_ *a, case_stmt *c, unsigned subcaseid, double factor);

#endif /* CASE_H */

