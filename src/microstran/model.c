#define MSTRANP_BUILD

#include "model.h"
#include "displacements.h"
#include "new.h"
#include "case.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

//#define NODISPLACEMENT_WARN

node_stmt node_create(unsigned id,double x,double y,double z,unsigned flags){
	node_stmt n;
	n.id = id;
	n.x = x;
	n.y = y;
	n.z = z;
	n.flags = flags;
	return n;
}

int node_print(FILE *f, const node_stmt *n){
	return fprintf(f,"NODE %5d %10f %10f %10f %06d\n",n->id, n->x, n->y, n->z, n->flags);
}

node_stmt *node_translate(node_stmt *n, double dx, double dy, double dz){
	n->x += dx;
	n->y += dy;
	n->z += dz;
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
	return m;
}

void model_destroy(model *a){
	
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

cbool model_add_prop(model *a, unsigned id, char libr[], char name[], char desc[]
		, cbool isdefault, double vals[MAXPROPVALS]
){
	prop_stmt p;
	//int i = 0;

	if(a->num_props > MAXPROPS)return 0;

	p.id = id; /* TODO check unique */
	strncpy(p.libr,libr,MAXPROPLIBNAME);
	strncpy(p.name,name,MAXPROPNAME);
	strncpy(p.desc,desc,MAXPROPDESC);
	p.isdefault = isdefault;

#ifdef MODEL_DEBUG
	for(i=0; i < MAXPROPVALS; ++i){
		p.vals[i] = vals[i];
	}
	fprintf(stderr,"PROP %d ",p.id);
	if(strlen(libr)){
		fprintf(stderr,"LIBR %s ",libr);
	}
	fprintf(stderr,"%s Y",name);
	if(isdefault){
		fprintf(stderr,"%s",(isdefault?"default":""));
	}
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
	case_stmt *c1;
	c1 = model_find_case(a,c->id);
	if(c1){
		fprintf(stderr,"case ID %d already in use\n",c->id);
		return 0;
	}

	array_append(&(a->cases), c);
	//fprintf(stderr,"Added case %d '%s' (%d %s)\n",c->id, c->name, c->data.num, (c->type==CASE_LOADS ? "loads" : "subcases"));
	return 1;
}

case_stmt *model_find_case(model *a,unsigned caseid){
	unsigned i,n;
	case_stmt *c;
	n = ARRAY_NUM(a->cases);
	for(i=0; i < n; ++i){
		c = array_get(&(a->cases), i);
		//fprintf(stderr,"looking at case %d, id = %d\n",i,c->id);
		if(c->id == caseid){
			return c;
		}
	}
	return NULL;
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

#define sq(X) ((X)*(X))
		l = sqrt(sq(to->x - from->x) + sq(to->y - from->y) + sq(to->z - from->z));
		
		printf("%d\t%d\t%d\t%d\t%s\t%f\n", memb->id, from->id, to->id, prop->id, prop->name, l);
	}
}



