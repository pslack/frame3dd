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

#include "model.h"
#include "vec3.h"
#include "displacements.h"
#include "new.h"
#include "case.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* #define NODISPLACEMENT_WARN */
#define MSTRANP_NORMALISE_ORIENTATION
/* #define CASE_DEBUG */

#define PI 3.14159265358

node_stmt node_create(unsigned id,vec3 pos,unsigned flags){
	node_stmt n;
	n.id = id;
	n.pos = pos;
	n.flags = flags;
	return n;
}

int node_print(FILE *f, const node_stmt *n){
#define BIT(VAL) (n->flags & MSTRANP_NODE_FIX##VAL ? '1' : '0')
	return fprintf(f,"NODE %5d %10f %10f %10f %c%c%c%c%c%c\n",n->id, n->pos.x, n->pos.y, n->pos.z
		, BIT(X), BIT(Y), BIT(Z)
		, BIT(MX), BIT(MY), BIT(MZ)
	);
#undef BIT
}

node_stmt *node_translate(node_stmt *n, double dx, double dy, double dz){
	n->pos.x += dx;
	n->pos.y += dy;
	n->pos.z += dz;
	return n;
}

cbool model_add_node(model *a, node_stmt node){
	//fprintf(stderr,"num nodes = %d\n",a->num_nodes);
	if(a->num_nodes > MAXNODES)return 0;
	/* TODO check ID unique */
	a->node[a->num_nodes++] = node;
	//fprintf(stderr,"NODE %d %f %f %f %d\n",node.id,node.x,node.y,node.z,node.flags);
	return 1;
}

model *model_create(unsigned version, unsigned type, unsigned vert, unit_stmt unit){
	model *a;
	a = NEW(model);
	a->version = version;
	a->type = type;
	a->vert = vert;
	a->unit = unit;

	a->num_nodes = 0;
	a->num_membs = 0;
	a->num_props = 0;
	a->num_matls = 0;
	a->cases = ARRAY_CREATE(case_stmt, 10);
	a->moffs = ARRAY_CREATE(moff_stmt, 10);
	return a;
}

model *model_copy(const model *om){
	model *m = NEW(model);
	if(!memcpy(m,om,sizeof(model))){
		fprintf(stderr,"unable to copy model\n");
		abort();
	}
	// case_stmts need special treatment due to dynamic memory allocation
	unsigned i;
	case_stmt *c;
	m->cases = ARRAY_CREATE(case_stmt,ARRAY_NUM(om->cases));
	for(i=0; i < ARRAY_NUM(om->cases); ++i){
		c = case_copy(array_get(((array *)&(om->cases)),i));
		array_set(&(m->cases), i, c);
		free(c); /* data at 'c' will have been copied now */
	}

	// moff_stmts need special treatment due to dynamic memory allocation
	moff_stmt moff;
	m->moffs = ARRAY_CREATE(moff_stmt,ARRAY_NUM(om->moffs));
	for(i=0; i < ARRAY_NUM(om->moffs); ++i){
		moff = *((moff_stmt *)array_get((array *)&(om->moffs),i));
		array_set(&(m->moffs), i, &moff);
	}

	return m;
}

void model_destroy(model *a){
	array_destroy(&(a->cases));
	array_destroy(&(a->moffs));
	free(a);
}

cbool model_find_node(const model *a, unsigned id, unsigned *index){
	int i;
	for(i=0; i <a->num_nodes; ++i){
		if(a->node[i].id == id){
			*index = i;
			return 1;
		}
	}
	return 0;
}

cbool model_add_memb(model *a, unsigned id,unsigned fromnode
		,unsigned tonode, member_orientation orient, unsigned prop, unsigned matl
		,unsigned flags1,unsigned flags2
){
	unsigned index;
	memb_stmt m;
	if(a->num_membs > MAXMEMBS)return 0;
	if(model_find_node(a,fromnode,&index)){
		m.fromnode = index;
	}else{
		fprintf(stderr,"Invalid MEMB 'from' number\n");
		return 0;
	}
	if(model_find_node(a,tonode,&index)){
		m.tonode = index;
	}else{
		fprintf(stderr,"Invalid MEMB 'to' number\n");
		return 0;
	}
	if(orient.axis=='\0'){
		if(model_find_node(a,orient.node, &index)){
			orient.node = index;
			/* fprintf(stderr,"Member %d is oriented to node %d\n",id,a->node[orient.node].id); */
		}else{
			fprintf(stderr,"Invalid MEMB orientation node '%d'\n",orient.node);
			return 0;
		}
	}
	m.id = id; /* TODO check unique */
	m.orient = orient;
	m.prop = prop;
	m.matl = matl;
	m.flags1 = flags1;
	m.flags2 = flags2;
	//fprintf(stderr,"MEMB %d %d %d %c%c %d %d %d %d\n",m.id,fromnode,tonode,axisdir,axis,prop,matl,flags1,flags2);
	a->memb[a->num_membs++] = m;
	return 1;
}

cbool model_find_memb(const model *a, const unsigned membid, unsigned *membindex){
	unsigned i;
	if(a->num_membs==0){
		fprintf(stderr,"No MEMBs in model!\n");
		return 0;
	}
	for(i=0; i<a->num_membs; ++i){
		if(a->memb[i].id == membid){
			*membindex = i;
			return 1;
		}
	}
	fprintf(stderr,"No member with ID %d was found.\n", membid);
	return 0;
}


cbool model_find_memb_from_to(const model *a, const unsigned fromnodeid, const unsigned tonodeid, unsigned *membindex){
	unsigned fromnodeindex, tonodeindex,i;
	if(a->num_membs==0){
		fprintf(stderr,"No MEMBs in model!\n");
		return 0;
	}
	if(!model_find_node(a,fromnodeid,&fromnodeindex)){
		fprintf(stderr,"Invalid 'fromnodeid' in model_find_memb_from_to\n");
		return 0;
	}
	if(!model_find_node(a,tonodeid,&tonodeindex)){
		fprintf(stderr,"Invalid 'tonodeid' in model_find_memb_from_to\n");
		return 0;
	}
	for(i=0; i<a->num_membs; ++i){
		if((a->memb[i].fromnode == fromnodeindex && a->memb[i].tonode == tonodeindex)
			|| (a->memb[i].tonode == fromnodeindex && a->memb[i].fromnode == tonodeindex)
		){
			*membindex = i;
			return 1;
		}
	}
	fprintf(stderr,"No MEMB between these two nodes\n");
	return 0;
}

MSTRANP_API memb_stmt *model_find_memb_from(model *a, const unsigned nodeid1, memb_stmt *start){
	memb_stmt *i = start;
	if(start==NULL){
		i = a->memb;
		i--;
		/* now start = a->memb[-1] */
	}
	memb_stmt *end = &a->memb[a->num_membs];
	while(++i < end){
		if(
			a->node[i->fromnode].id == nodeid1
			|| a->node[i->tonode].id == nodeid1
		){
			return i;
		}
	}
	return NULL;
}

/**
	Member orientation is defined by a third node ID, or
	else several axis orientations X, Y, Z or -X, -Y, -Z.

	@TODO FIXME still need to sort out 'D' orientation
	as specified for Microstran file format.
*/
vec3 memb_get_orientation(const model *a, const memb_stmt *m){
	vec3 AC = vec3_create(0,0,0);
	vec3 AB = vec3_diff(a->node[m->tonode].pos, a->node[m->fromnode].pos);
	switch(m->orient.axis){
		case 'X': AC.x = 1; break;
		case 'Y': AC.y = 1; break;
		case 'Z': AC.z = 1; break;
		case '\0':
			AC = vec3_diff(a->node[m->orient.node].pos, a->node[m->fromnode].pos);
			break;
		default:
			fprintf(stderr,"Invalid orient.axis value\n");
			exit(1);
	}

	if(m->orient.axis!='\0'){
		if(m->orient.dir=='-'){
			AC = vec3_scale(AC,-1);
		}
	}

#if 0
	fprintf(stderr,"axis = %c\n",m->orient.axis);
	fprintf(stderr,"\nAB = ");
	vec3_print(stderr,AB);
	fprintf(stderr,"\nAC = ");
	vec3_print(stderr,AC);
	fprintf(stderr,"\n");
#endif

#ifdef MSTRANP_NORMALISE_ORIENTATION
	return vec3_norm(vec3_cross(AC,AB));
#else
	return vec3_cross(AC,AB);
#endif

#if 0
	fprintf(stderr,"O = ");
	vec3_print(stderr,O);
	fprintf(stderr,"\n");
#endif

}

cbool model_add_member_offset(model *a, unsigned memberid
		, coord_sys_t coordsys, vec3 deltafrom, vec3 deltato
){
	/* MOFF statements are see later in the .arc file, and we're storing them
	in a different bit of memory for compatibility with the apparent philosophy
	of the Microstran approach. All the same, this is still a rather naive
	implementation and will hopefully be improved once we understand its
	workings a bit better. */
	moff_stmt m;
	m.id = memberid;
	m.coordsys = coordsys;
	m.deltafrom = deltafrom;
	m.deltato = deltato;
	array_append(&(a->moffs), &m);
	return 1; /* success */
}

moff_stmt *model_find_member_offset(const model *a, const unsigned memberid){
	unsigned i, n;
	n = ARRAY_NUM(a->moffs);
	for(i=0; i<n; ++i){
		moff_stmt *moff = (moff_stmt *)array_get((array *)&(a->moffs), i);
		if(moff->id == memberid){
			return moff;
		}
	}
	return NULL;
}

#if 0
/* this code replaced by use of ctrans_apply in arc2iv.cpp */
cbool model_get_member_offset_global(const model *a, const unsigned memberid, moff_stmt *moff){
	cbool found=0;
	moff_stmt o;
	const memb_stmt *m;
	vec3 X;
	unsigned i, n;
	const node_stmt *A, *B;

	assert(moff!=NULL);

	n = ARRAY_NUM(a->moffs);
	for(i=0; i<n; ++i){
		// copy the member offset to the
		o = *(moff_stmt *)array_get((array *)&(a->moffs), i);
		if(o.id == memberid){
			found=1;
			break;
		}
	}
	if(!found)return 0;

	moff_print(stderr, &o);

	if(o.coordsys==MSTRANP_COORDS_GLOBAL){
		*moff = o;
	}else if(o.coordsys==MSTRANP_COORDS_LOCAL){
		unsigned memberindex;
		if(!model_find_memb(a, o.id, &memberindex))return 0;
		m = &(a->memb[memberindex]);

		moff->id = o.id;
		moff->coordsys = MSTRANP_COORDS_GLOBAL;

		X = memb_get_orientation(a, m);
		fprintf(stderr,"mod(X) = %f\n", vec3_mod(X));

		A = &(a->node[m->fromnode]);
		B = &(a->node[m->tonode]);

		vec3 AB = vec3_norm(vec3_diff(B->pos, A->pos));
		vec3 zzdir;
		double theta_zz = vec3_angle_cross(vec3_create(0,0,1), AB, &zzdir);
		if(vec3_mod(zzdir)<1e-6){
			zzdir = vec3_create(0,0,1);
			theta_zz = 0;
		}

		assert(!isnan(theta_zz));
		assert(!isnan(zzdir.x));

		vec3 Xdd = vec3_rotate(X,zzdir, -theta_zz);
		assert(fabs(vec3_mod(Xdd)-1.0)<1e-6);

		vec3 xxxdir;
		double theta_xxx = vec3_angle_cross(vec3_create(1,0,0), Xdd, &xxxdir);

		if(vec3_mod(xxxdir)<1e-6){
			fprintf(stderr,"xxxdir = 0\n");
			xxxdir = vec3_create(0,0,1);
			theta_xxx = 0;
		}

		assert(!isnan(theta_xxx));
		assert(!isnan(xxxdir.x));

		fprintf(stderr,"LO deltafrom = %f %f %f\n",o.deltafrom.x, o.deltafrom.y, o.deltafrom.z);
		fprintf(stderr,"rotate by %f deg about %f %f %f\n", -theta_xxx*180./PI, xxxdir.x, xxxdir.y, xxxdir.z);

		vec3 Add = vec3_rotate(o.deltafrom,xxxdir,-theta_xxx);
		assert(!isnan(Add.x));
		assert(!isnan(Add.y));
		assert(!isnan(Add.z));

		/* we can now perform the coordinate transform via rotation xxx then zz */
		moff->deltafrom = vec3_rotate(Add,zzdir,-theta_zz);

		vec3 Bdd = vec3_rotate(o.deltato,xxxdir,-theta_xxx);
		moff->deltato = vec3_rotate(Bdd,zzdir,-theta_zz);

		moff_print(stderr, moff);

		assert(!isnan(moff->deltafrom.x));
		assert(!isnan(moff->deltafrom.y));
		assert(!isnan(moff->deltafrom.z));
		assert(!isnan(moff->deltato.x));
		assert(!isnan(moff->deltato.y));
		assert(!isnan(moff->deltato.z));

	}
	return 1;
}
#endif

int moff_print(FILE *f, const moff_stmt *o){
	const char *code = "??";
	switch(o->coordsys){
		case MSTRANP_COORDS_GLOBAL: code="GL"; break;
		case MSTRANP_COORDS_LOCAL:  code="LO"; break;
	}
	return fprintf(f,"MOFF %5d %s %10f %10f %10f %10f %10f %10f\n",o->id, code
		, o->deltafrom.x, o->deltafrom.y, o->deltafrom.z
		, o->deltato.x, o->deltato.y, o->deltato.z
	);
}


cbool model_add_prop(model *a, unsigned id, char libr[], char name[], char desc[]
		, double vals[MAXPROPVALS]
){
	prop_stmt p;
	//int i = 0;

	if(a->num_props > MAXPROPS)return 0;

	p.id = id; /* TODO check unique */
	strncpy(p.libr,libr,MAXPROPLIBNAME);
	strncpy(p.name,name,MAXPROPNAME);
	strncpy(p.desc,desc,MAXPROPDESC);

#ifdef MODEL_DEBUG
	for(i=0; i < MAXPROPVALS; ++i){
		p.vals[i] = vals[i];
	}
	fprintf(stderr,"PROP %d ",p.id);
	if(strlen(libr)){
		fprintf(stderr,"LIBR %s ",libr);
	}
	fprintf(stderr,"%s Y",name);
	fprintf(stderr,"\n\t");
	for(i=0; i<MAXPROPVALS; ++i){
		fprintf(stderr,"%e ",vals[i]);
	}
	fprintf(stderr,"\n");
#endif
	a->prop[a->num_props++] = p;
	return 1;
}

cbool model_add_matl(model *a, unsigned id, double E, double sigma_y
		, double rho, double beta
){
	matl_stmt m;
	if(a->num_matls > MAXMATLS)return 0;
	m.id = id; /* TODO check unique */
	m.E = E;
	m.sigma_y = sigma_y;
	m.rho = rho;
	m.beta = beta;
	//fprintf(stderr,"MATL %d %e %e %e %e\n",m.id,m.E,m.sigma_y,m.rho,m.beta);
	a->matl[a->num_matls++] = m;
	return 1;
}

cbool model_add_case(model *a,case_stmt *c){
	unsigned caseindex;
	if(model_find_case(a,c->id, &caseindex)){
		fprintf(stderr,"case ID %d already in use\n",c->id);
		return 0;
	}

	array_append(&(a->cases), c);
	if(!model_find_case(a,c->id, &caseindex)){
		fprintf(stderr,"Failed to add case %d\n",c->id);
		return 0;
	}

#ifdef CASE_DEBUG
	/* case_stmt *c1; */
	c1 = array_get(&a->cases, caseindex);
	fprintf(stderr,"Added case %d '%s' (%d %s)\n",c1->id, c1->name, ARRAY_NUM(c1->data), (c1->type==CASE_LOADS ? "loads" : "subcases"));
	if(c1->type ==CASE_COMB){
		unsigned i;
		for(i=0; i<ARRAY_NUM(c1->data); ++i){
			comb_stmt *co = array_get(&c1->data, i);
			case_stmt *c2 = array_get(&a->cases, co->caseindex);
			fprintf(stderr," %c %f × CASE %d '%s'\n",(i==0 ? '=' : '+'), co->factor, c2->id, c2->name);
		}
	}
#endif
	return 1;
}

cbool model_find_case(model *a,unsigned caseid, unsigned *caseindex){
	unsigned i,n;
	case_stmt *c;
	n = ARRAY_NUM(a->cases);
	for(i=0; i < n; ++i){
		c = array_get(&(a->cases), i);
		//fprintf(stderr,"looking at case %d, id = %d\n",i,c->id);
		if(c->id == caseid){
			*caseindex = i;
			return 1;
		}
	}
	return 0;
}

case_stmt *model_get_case(model *a, unsigned caseindex){
	if(caseindex >= ARRAY_NUM(a->cases))return NULL;
	return (case_stmt *)array_get(&a->cases,caseindex);
}

cbool model_apply_displacements(model *m, casedisplacements *cd){
	/* important note: nodes are not densely numbered, so we must use 'find' functions etc */
	/* note that we required all nodes to be given a displacement */
	unsigned im, nodeid;
	nodedisplacement *nd;
	unsigned displaced = 0;

	double maxdisp = 0, disp;
	unsigned nodemaxdisp = 999999999;

	fprintf(stderr,"Displacing model with %d nodes...\n",m->num_nodes);
	//fprintf(stderr,"Displacement set has %d nodes...\n",ARRAY_NUM(cd->nodes));

	for(im=0; im < m->num_nodes; ++im){
		nodeid = m->node[im].id;
		//fprintf(stderr,"Model node %d\n",nodeid);
		/* find the required displacement for this node */
		nd = casedisplacements_find_node(cd, nodeid);
		if(!nd){
#ifdef NODISPLACEMENT_WARN
			fprintf(stderr,"Warning: displacement for node %d not found\n",nodeid);
#endif
			continue;
		}
		node_translate(&(m->node[im]), nd->dx, nd->dy, nd->dz);
		disp = sqrt(nd->dx*nd->dx + nd->dy*nd->dy + nd->dz*nd->dz);
		if(disp>maxdisp){
			maxdisp=disp;
			nodemaxdisp=nodeid;
		}
		displaced++;
	}

	fprintf(stderr,"Displaced %d nodes (max displacement %f at node %d)\n",displaced,maxdisp,nodemaxdisp);
	return 1;
}

prop_stmt *model_find_prop(model *m, unsigned propid){
	unsigned i;
	for(i=0; i< m->num_props; ++i){
		if(m->prop[i].id == propid){
			return m->prop + i;
		}
	}
	return NULL;
}

void model_write_inventory(model *m){
	unsigned i;
	double l;
	memb_stmt *memb;
	prop_stmt *prop;
	node_stmt *from, *to;

	printf("MEMBER\tFROM\tTO\tPROPID\tSECTION\tLENGTH\n");
	for(i=0; i < m->num_membs; ++i){
		memb = m->memb + i;
		prop = model_find_prop(m, memb->prop);
		if(prop==NULL){
			fprintf(stderr,"UNKNOWN PROP ID");
			exit(1);
		}
		from = m->node + memb->fromnode;
		to = m->node + memb->tonode;

		l = vec3_mod(vec3_diff(to->pos,from->pos));

		printf("%d\t%d\t%d\t%d\t%s\t%f\n", memb->id, from->id, to->id, prop->id, prop->name, l);
	}
}



