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

#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>

#include <unistd.h>
#include <stdlib.h>

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

ostream & operator<< (ostream &os, const SbVec3f &v){
	os << v[0] << " " << v[1] << " " << v[2] << " ";
	return os;
}


double defaultscaleF = 5;
double defaultscaleM = 1;
const char *defaultmodelfname = "model.arc";
const char *defaultsceneoutfile = "forcebalance.iv";

void usage(char *progname){
	cerr << progname << ": Parse and display force balance on nodes using Microstran output." << endl;
	cerr << "Usage: " << progname << " [-MF] [-o[OUTFILE]] [-n NODEID] [-m MODEL.arc] [-f FORCES.p1]" << endl;
	cerr << "    -m MODEL    Filename of Microstran model .arc file to read. Defaults to '" << defaultmodelfname << "'." << endl;
	cerr << "    -f FORCES   Filename of .p1 member forces report to read. Defaults to MODEL.p1 (by changing file extension)." << endl;
	cerr << "    -s SCALEF   Scale forces in visual output, such that SCALEF kN = 1 m" << endl;
	cerr << "    -t SCALEM   Scale moments in visual output, such that SCALEM kNm = 1 m" << endl;
	cerr << "    -n NODEID   Show only forces for node NODEID" << endl;
	cerr << "    -o[OUTFILE] Filename to write to, in Open Inventor .iv format (instead of rendering to screen). If name omitted, defaults to '" << defaultsceneoutfile << "'." << endl;
	cerr << "    -M          Toggle display of moments (default off)" << endl;
	cerr << "    -F          Toggle display of forces (default on)" << endl;
	cerr << "    -X          Toggle display of axes (default off)" << endl;
	exit(1);
}

int main(int argc, char **argv){

	const char *modelfname = defaultmodelfname;

	const char *forcefile = NULL;
#define FILENMAX 200
	char forcefileauto[FILENMAX];

	const char *sceneoutfile = NULL;

	unsigned caseid = 1;

	bool showall = 1;
	unsigned shownode = 0;

	double scaleF = defaultscaleF;
	double scaleM = defaultscaleM;

	bool showmoments = false;
	bool showforces = true;
	bool showaxes = false;
	bool showoffsets = false;
	bool applyoffsets = true;

	char c;
    while((c = getopt (argc, argv, "MFn:m:f:c:s:t:o::")) != -1){
		switch(c){
			case 'F': showforces = !showforces; break;
			case 'M': showmoments = !showmoments; break;
			case 'X': showaxes = !showaxes; break;
			case 'n': showall = false; shownode = atoi(optarg); break;
			case 'm': modelfname = optarg; break;
			case 'f': forcefile = optarg; break;
			case 'c': caseid = atoi(optarg); break;
			case 's': scaleF = atof(optarg); break;
			case 't': scaleM = atof(optarg); break;
			case 'o':
				if(optarg)sceneoutfile = optarg;
				else sceneoutfile = defaultsceneoutfile;
				break;
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

	if(forcefile == NULL){
		strncpy(forcefileauto,(char *)modelfname,FILENMAX);
		char *c = forcefileauto;
		while(*c!='\0')++c;
		while(*c!='.')--c;
		if(c==modelfname)throw runtime_error("Unable to generate forces filename, please specify using '-f' option");
		sprintf(c,".p1");
		forcefile = forcefileauto;
	}


	cerr << "Loading member forces file '" << forcefile << "'..." << endl;
	modelforces *MF = NULL;
	f1 = fopen(forcefile,"r");
	if(!f1){
		cerr << "ERROR: Unable to read member forces file" << endl;
		exit(1);
	}

	fclose(f1);
	p = parseCreateFileName(forcefile);
	if(!parseMicrostranForces(p,&MF)){
		cerr << "ERROR: Unable to parse member forces file '" << forcefile << "'" << endl;
		exit(1);
	}

	fprintf(stderr,"\nParse member forces OK, %d load cases found\n",ARRAY_NUM(MF->cases));

	//---------------------
	// find the results for the load case we seek

	cerr << "Processing member forces for case " << caseid << "..." << endl;

	unsigned caseindex;
	if(!model_find_case(M,caseid,&caseindex)){
		cerr << "Unable to find load-case in model";
		exit(1);
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
    if(showaxes){
		root->addChild(axes(0.2));
	}

	double maxF = 0, maxM = 0;
	bool foundloads = false;
	bool foundfix = false;
	bool founddelta = false;
	double F_bal_max = 0; //< maximum force residual (should be close to zero)
	double M_bal_max = 0; //< maximum moment residual (should be close to zero)
	unsigned node_F_bal_max = 0;
	unsigned node_M_bal_max = 0;

	if(!applyoffsets){
		cerr << "NOTE: moments due to member offsets are NOT being applied" << endl;
	}

	// loop through all the nodes
	for(unsigned i = 0; i < M->num_nodes; ++i){
		unsigned nodeid = M->node[i].id;

		ndld_stmt nl;
		if(case_total_load_node(M, cl, nodeid, &nl)){
			foundloads = true;
			if(showall || nodeid==shownode){
				if(showforces){
					fprintf(stderr, "Node %d: ",nodeid);
					VEC3_PR(nl.F);
				}
				if(showmoments){
					fprintf(stderr, "Node %d: ",nodeid);
					VEC3_PR(nl.M);
				}
			}
		}

		unsigned nodeindex;
		if(!model_find_node(M, nodeid, &nodeindex)){
			throw runtime_error("Unknown node ID");
		}
		vec3 pos = M->node[nodeindex].pos;

		if(showall || nodeid==shownode){
			if(showforces){
				stringstream ss;
				ss << "F=" << vec3_mod(nl.F);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(nl.F,1./scaleF)),GREEN,ss.str().c_str(),0.01));
			}
			if(showmoments){
				stringstream ss;
				ss << "M=" << vec3_mod(nl.M);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(nl.M,1./scaleM)),PURPLE,ss.str().c_str(),0.01));
			}
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
			if(showaxes &&
					(showall || thisnode->id==shownode)
			){
				stringstream ss;
				ss << "X-ax " << m->id << ": " << fromnode->id << " -- " << tonode->id;
				root->addChild(arrow(
					from_vec3(fromnode->pos)
					,from_vec3(vec3_add(fromnode->pos, vec3_scale(vec3_norm(minax),0.3)))
					,LIME
					,ss.str().c_str()
					,0.005
				));
			}

			ctrans_matrix ct = ctrans_rotation_axes(vec3_diff(tonode->pos,fromnode->pos),minax);

			vec3 M_gl = vec3_negate(ctrans_apply(ct,M_local));
			vec3 F_gl = ctrans_apply(ct,F_local);

			if(showmoments && !showall && thisnode->id==shownode){
				cerr << "Moment due to member " << m->id << " (from node " << othernode->id << "): " << from_vec3(M_gl) << endl;
			}

			if(vec3_mod(F_gl)>maxF)maxF = vec3_mod(F_gl);
			if(vec3_mod(M_gl)>maxM)maxM = vec3_mod(M_gl);

			if((showall || thisnode->id == shownode) && vec3_mod(F_gl) > 5e-4){
				if(showforces){
					stringstream ss;
					ss << othernode->id << ":F ";
					ss << vec3_mod(F_gl);
					root->addChild(arrow(
						from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(F_gl,1./scaleF))
						,(m->id >= 7200 ? YELLOW : CYAN)
						,ss.str().c_str()
						,0.01
					));
				}
				if(showmoments){
					stringstream ss;
					ss << othernode->id << ": ";
					ss << vec3_mod(M_gl);
					root->addChild(arrow(
						from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(M_gl,1./scaleM))
						,(m->id >= 7200 ? MAGENTA : PURPLE)
						,ss.str().c_str()
						,0.01
					));
				}
			}

			M_res = vec3_add(M_res, M_gl);
			F_res = vec3_add(F_res, F_gl);

			// add offset forces to resultant moment
			if(delta){
				founddelta = true;
				vec3 delta_gl = ctrans_apply(ct, *delta);
				vec3 M_offset = vec3_cross(delta_gl, F_gl);
				if(showoffsets){
					if(showall || nodeid == shownode){
						cerr << "Offset for member " << m->id << " = ";
						vec3_print(stderr,delta_gl);
						cerr << " --> M_offset = " << vec3_mod(M_offset) << " = |";
						vec3_print(stderr,M_offset);
						cerr << "|" << endl;
					}else if(nodeid == shownode){
						cerr << "No offset for member " << m->id << endl;
					}
				}
				if(applyoffsets){
					M_res = vec3_add(M_res, M_offset);
				}
			}
		}

		vec3 F_bal = vec3_add(F_res, nl.F);
		vec3 M_bal = vec3_add(M_res, nl.M);

		if(showall || nodeid == shownode){
			cerr << "Resultant force at node " << nodeid << " = ";
			vec3_print(stderr,F_res);
			cerr << endl;

#ifdef SHOW_TOTAL_REACTION
			if(vec3_mod(F_res) > 5e-4){
				stringstream ss;
				ss << "TR:";
				ss << vec3_mod(F_res);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(F_res,1./scaleF)),YELLOW,ss.str().c_str()));
			}
#endif

			if(showforces){
				cerr << "Force error at node " << nodeid << " = " << vec3_mod(F_bal) << " = |";
				vec3_print(stderr,F_bal);
				cerr << "|" << endl;

				stringstream ss;
				ss << "B:";
				ss << vec3_mod(F_bal);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(F_bal,1./scaleF)),PURPLE,ss.str().c_str(),0.01));
			}

			if(showmoments){
				cerr << "Moment error at node " << nodeid << " = " << vec3_mod(M_bal) << " = |";
				vec3_print(stderr,M_bal);
				cerr << "|" << endl;

				stringstream ss;
				ss << "M_bal:";
				ss << vec3_mod(M_bal);
				root->addChild(arrow(from_vec3(pos), from_vec3(pos) + from_vec3(vec3_scale(M_bal,1./scaleM)),ORANGE,ss.str().c_str(),0.01));
			}
		}

		// to calculate the maximum force/moment errors, cancel out restrained components

		bool cleared = false;
		if(M->node[i].flags & MSTRANP_NODE_FIXX){
			if(showforces && (showall || nodeid == shownode))cerr << "Node " << nodeid << " has 'X' restraint" << endl;
			cleared = true; foundfix = true;
			F_bal.x = 0;
		}
		if(M->node[i].flags & MSTRANP_NODE_FIXY){
			if(showforces&& (showall || nodeid == shownode))cerr << "Node " << nodeid << " has 'Y' restraint" << endl;
			cleared = true; foundfix = true;
			F_bal.y = 0;
		}
		if(M->node[i].flags & MSTRANP_NODE_FIXZ){
			if(showforces && (showall || nodeid == shownode))cerr << "Node " << nodeid << " has 'Z' restraint" << endl;
			cleared = true; foundfix = true;
			F_bal.z = 0;
		}
		if(M->node[i].flags & MSTRANP_NODE_FIXMX){
			if(showmoments && (showall || nodeid == shownode))cerr << "Node " << nodeid << " has 'MX' restraint" << endl;
			cleared = true;
			M_bal.x = 0;
		}
		if(M->node[i].flags & MSTRANP_NODE_FIXMY){
			if(showmoments && (showall || nodeid == shownode))cerr << "Node " << nodeid << " has 'MY' restraint" << endl;
			cleared = true;
			M_bal.y = 0;
		}
		if(M->node[i].flags & MSTRANP_NODE_FIXMZ){
			if(showmoments && (showall || nodeid == shownode))cerr << "Node " << nodeid << " has 'MZ' restraint" << endl;
			cleared = true;
			M_bal.z = 0;
		}
		if(cleared && (showall || nodeid == shownode)){
			node_print(stderr,&M->node[i]);
		}

		if(vec3_mod(F_bal)>F_bal_max){
			F_bal_max = vec3_mod(F_bal);
			node_F_bal_max = nodeid;
		}
		if(vec3_mod(M_bal)>M_bal_max){
			M_bal_max = vec3_mod(M_bal);
			node_M_bal_max = nodeid;
		}
	}

	if(!foundloads){
		cerr << "WARNING: no applied node-loads found!" << endl;
	}
	if(!foundfix){
		cerr << "WARNING: no translational restraints found!" << endl;
	}
	if(!founddelta){
		cerr << "NOTE: no member offsets found!" << endl;
	}
	cerr << "Largest member force = " << maxF << endl;
	cerr << "Largest member moment = " << maxM << endl;
	cerr << "Maximum force error = " << F_bal_max << " (at node " << node_F_bal_max << ")" << endl;
	cerr << "Maximum moment error = " << M_bal_max << " (at node " << node_M_bal_max << ")" << endl;

	if(sceneoutfile){
		SoWriteAction wr;
		FILE *f;
		f = fopen(sceneoutfile,"w");
		if(!f) throw runtime_error("Unable to write to scene output file");
		wr.getOutput()->setFilePointer(f);
		cerr << "Writing scene file to '" << sceneoutfile << "'...";
		wr.apply(root);
		cerr << "done!" << endl;
	}else{
		fprintf(stderr,"Rendering scene to OpenGL directly...\n");

		// Use one of the convenient SoQt viewer classes.
		SoQtExaminerViewer *eviewer = new SoQtExaminerViewer(mainwin);
		eviewer->getGLRenderAction()->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_BLEND);
		eviewer->setSceneGraph(root);
		eviewer->show();

		// Pop up the main window.
		SoQt::show(mainwin);
		// Loop until exit.
		SoQt::mainLoop();

		// Clean up resources.
		delete eviewer;
	}
	root->unref();

	exit(0);
}

