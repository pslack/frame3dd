/** @FILE
	This program displays the forces applied at each node in the structure
	and creates a graphical output of the forces on each node from each
	of the node's connected members. The total reaction force on the
	node is shown, as well as the total applied load (for the selected load
	case), and finally the force balance (which should be very close to zero).

	You can set the program to display forces for all nodes in the structure,
	or your can specify a single node, using a command-line argument.
*/

#include <microstran/config.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>

#include <microstran/parse.h>
#include <microstran/model.h>
#include <microstran/modelparser.h>
#include <microstran/forceparser.h>
#include <microstran/ctrans.h>
#undef assign
#undef done


#include "render.h"
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/actions/SoGLRenderAction.h>

using namespace std;

static SbVec3f from_vec3(vec3 A){
	return SbVec3f(A.x, A.y, A.z);
}

#if 0
Vector::operator()(const vec3 &v){
	return Vector(v.x,v.y,v.z);
}
#endif

double defaultscaleF = 5;

void usage(char *progname){
	cerr << progname << ": process front-surface join forces." << endl;
	cerr << "Usage: " << progname << " [-n NODEID] [-m MODEL.arc] [-f FORCES.p1]" << endl;
	cerr << "    -m FILE   Microstran model .arc file to read." << endl;
	cerr << "    -f FILE   Member forces data file to read (.p1 report format from Microstran)" << endl;
	cerr << "    -s SCALE  Scale forces in visual output, such that SCALE kN = 1 m" << endl;
	cerr << "    -n NODEID Show only forces for node NODEID" << endl;
	exit(1);
}

int main(int argc, char **argv){

	const char *defaultmodelfname = "model.arc";
	const char *modelfname = defaultmodelfname;

	const char *defaultforcefile = "forces.p1";
	const char *forcefile = defaultforcefile;

	unsigned caseid = 1;

	bool showall = 1;
	unsigned shownode = 0;

	double scaleF = defaultscaleF, scaleM = 5;

	char c;
    while((c = getopt (argc, argv, "n:m:f:c:s:")) != -1){
		switch(c){
			case 'n': showall = false; shownode = atoi(optarg); break;
			case 'm': modelfname = optarg; break;
			case 'f': forcefile = optarg; break;
			case 'c': caseid = atoi(optarg); break;
			case 's': scaleF = atof(optarg); break;
			default:
				usage(argv[0]);
		}
	}

	//--------------------
	// read in microstran model

	parse *p;
	model *M;

	cerr << "Loading Microstran model '" << modelfname << "'" << endl;
	FILE *f1 = fopen(modelfname,"r");
	if(!f1)throw runtime_error("Unable to read Microstran model file");
	fclose(f1);
	p = parseCreateFileName(modelfname);
	if(parseModelMicrostran(p,&M)){
		fprintf(stderr,"\nParsed OK\n%d nodes, %d members, %d section properties, %d materials read.\n"
			,M->num_nodes,M->num_membs,M->num_props,M->num_matls
		);
	}else{
		fprintf(stderr,"\n\nPARSE FAILED\n");
		return 1;
	}
	parseDispose(p);

	//--------------------
	// read in member forces

	cerr << "Loading member forces file '" << forcefile << "'..." << endl;
	modelforces *MF = NULL;
	f1 = fopen(forcefile,"r");
	if(!f1)throw runtime_error("Unable to read member forces file");
	fclose(f1);
	p = parseCreateFileName(forcefile);
	if(!parseMicrostranForces(p,&MF))throw runtime_error("Unable to parse member forces file");
	fprintf(stderr,"\nParse member forces OK, %d load cases found\n",ARRAY_NUM(MF->cases));

	//---------------------
	// find the results for the load case we seek

	cerr << "Processing member forces for case " << caseid << "..." << endl;

	unsigned caseindex;
	if(!model_find_case(M,caseid,&caseindex)){
		throw runtime_error("Unable to find load-case in model");
	}
	case_stmt *cl = (case_stmt *)array_get(&M->cases, caseindex);

	cerr << "Case " << caseid  << ": '" << cl->name << "'" << endl;

#if 1
	caseforces *cf = modelforces_find_case(MF, caseid);
	if(!cf){
		throw runtime_error("Unable to find load-case in member forces data");
	}
#endif

	//----------------------
	// set up the scene for rendering forces graphically

    QWidget *mainwin = SoQt::init(argc, argv, argv[0]);

    SoSeparator *root = new SoSeparator;
    root->ref();
	root->addChild(axes());

	double maxF = 0, maxM = 0;
	bool foundloads = false;
	
	// loop through all the nodes
	for(unsigned i = 0; i < M->num_nodes; ++i){
		unsigned nodeid = M->node[i].id;

		ndld_stmt nl;
		if(case_total_load_node(M, cl, nodeid, &nl)){
			foundloads = true;
			fprintf(stderr, "Node %d: ",nodeid);
			VEC3_PR(nl.F);
		}

		unsigned nodeindex;
		if(!model_find_node(M, nodeid, &nodeindex)){
			throw runtime_error("Unknown node ID");
		}
		vec3 pos = M->node[nodeindex].pos;		

		if(1 || showall || nodeid==shownode){
			stringstream ss;
			ss << vec3_mod(nl.F);
			root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(nl.F,1./scaleF)),GREEN,ss.str().c_str()));
		}

		vec3 F_res = VEC3_ZERO;
		vec3 M_res = VEC3_ZERO;

		// iterate through members that are connected to the current one
		memb_stmt *m = NULL;
		while((m = model_find_memb_from(M, nodeid, m)) != NULL){
			//if(!(m->id >= 7200))continue;

			//cerr << "member = " << m->id << ": from " << M->node[m->fromnode].id << " to " << M->node[m->tonode].id << endl;

//			cerr << " + member" << endl;
			node_stmt *fromnode = &M->node[m->fromnode];
			node_stmt *tonode = &M->node[m->tonode];
			node_stmt *thisnode = NULL;
			node_stmt *othernode = NULL;

			memberforce *mf = caseforces_find_member(cf, m->id);

			assert(mf->member == m->id);		

			// check for member offset
			moff_stmt *moff = model_find_member_offset(M, m->id);
			vec3 *delta = NULL;

			membernodeforce *mnf = NULL;

			// force and mom of the member on the node, in member local coors
			vec3 F_local, M_local;

			assert(mf->tonode.node == tonode->id);
			assert(mf->fromnode.node = fromnode->id);

			if(fromnode->id == nodeid){
				thisnode = &M->node[m->fromnode];
				othernode = &M->node[m->tonode];
				mnf = &mf->fromnode;
				if(moff)delta = &moff->deltafrom;
				/*
					on the 'from' node, the force is negated twice: once
					for the fact of the sign convention, another for the
					fact that we are concerned with the force of the 
					member on the node, not vice versa
				*/
				F_local = mnf->F;
				M_local = mnf->M;
			}else if(tonode->id == nodeid){
				thisnode = &M->node[m->tonode];
				othernode = &M->node[m->fromnode];
				mnf = &mf->tonode;
				if(moff)delta = &moff->deltato;
				/* 
					on the 'to' node, the force is only negated once, for
					the fact of calculating force of member on node.
				*/
				F_local = vec3_negate(mnf->F);
				M_local = vec3_negate(mnf->M);
			}else{
				throw runtime_error("memberforce struct does not contain expected node (nodeid)");
			}

			//cerr << "'this' node is " << thisnode->id << endl;
			//VEC3_PR(F_local);
			//VEC3_PR(M_local);

			// get member orientation so that we can do the coordinate transformation
			vec3 minax = memb_get_orientation(M, m);
			if(showall || thisnode->id==shownode){
				stringstream ss;
				ss << "X-ax " << m->id << ": " << fromnode->id << " -- " << tonode->id;
				root->addChild(arrow(
					from_vec3(fromnode->pos)
					,from_vec3(vec3_add(fromnode->pos, vec3_scale(vec3_norm(minax),0.3)))
					,LIME
					,ss.str().c_str()
					,0.02
				));
			}

			ctrans_matrix ct = ctrans_rotation_axes(vec3_diff(tonode->pos,fromnode->pos),minax);

			vec3 M_gl = vec3_negate(ctrans_apply(ct,M_local));
			vec3 F_gl = ctrans_apply(ct,F_local);

			if(vec3_mod(F_gl)>maxF)maxF = vec3_mod(F_gl);
			if(vec3_mod(M_gl)>maxM)maxM = vec3_mod(M_gl);

			if((showall || thisnode->id == shownode) && vec3_mod(F_gl) > 5e-4){
				stringstream ss;
				ss << othernode->id << ": ";
				ss << vec3_mod(F_gl);
				root->addChild(arrow(
					from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(F_gl,1./scaleF))
					,(m->id >= 7200 ? YELLOW : CYAN)
					,ss.str().c_str()));
			}

			M_res = vec3_add(M_res, M_gl);
			F_res = vec3_add(F_res, F_gl);

			// add offset forces to resultant moment
			if(delta){
				vec3 delta_gl = ctrans_apply(ct, *delta);
				vec3 M_offset = vec3_cross(F_gl, delta_gl);
				M_res = vec3_add(M_res, M_offset);
			}
		}

		cerr << "Resultant force at node " << nodeid << " = ";
		vec3_print(stderr,F_res);
		cerr << endl;

		if(showall || nodeid == shownode){
#ifdef SHOW_TOTAL_REACTION
			if(vec3_mod(F_res) > 5e-4){
				stringstream ss;
				ss << "TR:";
				ss << vec3_mod(F_res);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(F_res,1./scaleF)),YELLOW,ss.str().c_str()));
			}
#endif
			{
				F_res = vec3_add(F_res, nl.F);
				stringstream ss;
				ss << "B:";
				ss << vec3_mod(F_res);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(F_res,1./scaleF)),PURPLE,ss.str().c_str()));
			}

#ifdef SHOW_MOMENTS
		if(vec3_mod(M_res) > 1e-3){
			root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(M_res,1./scaleM)),PURPLE));
		}
#endif
		}
	}

	if(!foundloads){
		cerr << "WARNING: no applied node-loads found!" << endl;
	}
	cerr << "Maximum force = " << maxF << endl;
	cerr << "Maximum moment = " << maxM << endl;

	SoWriteAction wr;
	const char *sceneoutfile = "forcebalance.iv";
	FILE *f;
	f = fopen(sceneoutfile,"w");
	if(!f) throw runtime_error("Unable to write to scene output file");
	wr.getOutput()->setFilePointer(f);
	cerr << "Writing scene file to '" << sceneoutfile << "'...";
	wr.apply(root);
	cerr << "done!" << endl;	

	exit(0);
}

