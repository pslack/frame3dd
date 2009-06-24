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
*//**@FILE
	Load and parse a microstran file, then render the node and member
	geometry using Coin3d / OpenGL.
*/

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/SoInput.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/actions/SoGLRenderAction.h>

#include <unistd.h>
#include <stdlib.h>
#include <limits.h>

#include <sstream>
#include <set>
#include <string>
#include <stdexcept>
#include <iostream>
using namespace std;

#include "render.h"

#include <microstran/modelparser.h>
#include <microstran/sectionsparser.h>
#include <microstran/ctrans.h>
#include <microstran/defaultpaths.h>

const char *defaultsceneoutfile = "microstranmodel.iv";
char defaultlibfile[1024];

void usage(const char *progname){
	fprintf(stderr,"Usage: %s [-o[OUTFILE]] [-a MIN] [-b MAX] [-m] [-h] [-l LIBFILE] [-t] INFILE\n",progname);
	fprintf(stderr,"Load and parse a Microstran .arc file and render using Coin3D/OpenGL\n");
	fprintf(stderr,"  -o[OUTFILE]  Open Inventor file to output (defaults to '%s'). No space after -o!\n",defaultsceneoutfile);
	fprintf(stderr,"  -a MEMBID    Ignore members with member IDs less than MEMBID.\n");
	fprintf(stderr,"  -b MEMBID    Ignore members with member IDs greater than MEMBID.\n");
	fprintf(stderr,"  -m           Ignore member offsets (don't apply the offers)\n");
	fprintf(stderr,"  -h           Render using higher quality graphics (slower).\n");
	fprintf(stderr,"  -l LIBFILE   Load a section library from LIBFILE.\n");
	fprintf(stderr,"  -t           Include text for node/member IDs and member sizes.\n");
	fprintf(stderr,"  INFILE       Microstran .arc file to render (eg 'model.arc')\n\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"If you specify -o, you will receive and output .iv file instead of an on-screen rendering.\n");
}

static SbVec3f vec3_to_coin(vec3 A){
	return SbVec3f(A.x,A.y,A.z);
}

int main(int argc, char **argv){

	const char *sceneoutfile = NULL;
	const char *libfile = defaultlibfile;
	bool highquality = false;
	bool infotext = false;
	bool memberoffsets = true;

	unsigned minmemb = 0;
	unsigned maxmemb = UINT_MAX;

	stringstream ss0;
	ss0 << get_default_data_path() << FRAME3DD_PATHSEP << "properties.txt";	
	strcpy(defaultlibfile,ss0.str().c_str());

	char c;
	while((c=getopt(argc,argv,"a:b:mho::l:t"))!=-1){
		switch(c){
			case 'a':
				minmemb = atoi(optarg);
				break;
			case 'b':
				maxmemb = atoi(optarg);
				break;
			case 'h':
				highquality = 1;
				break;
			case 'm':
				memberoffsets = 0;
				break;
			case 'o':
				sceneoutfile = (optarg ? optarg : defaultsceneoutfile);
				break;
			case 'l':
				libfile = optarg;
				break;
			case 't':
				infotext = true;
				break;
			case '?':
				usage(argv[0]);
				exit(1);
		}
	}

	if(optind != argc-1){
		fprintf(stderr,"%s: missing command-line argument (need a filename)\n",argv[0]);
		usage(argv[0]);
		exit(2);
	}
	const char *filename = argv[optind];

	parse *p1;
	fprintf(stderr,"Reading section library '%s'...",libfile);

	FILE *f;
	f = fopen(libfile,"r");
	if(f == NULL){
		fprintf(stderr,"\nUnable to open section library file '%s'!\n",libfile);
		exit(2);
	}
	fclose(f);

	p1 = parseCreateFileName(libfile);
	section_library *l = NULL;
	l = section_library_create();
	parseSections(p1,l);
	fprintf(stderr,"done!\n");
	parseDispose(p1);

	model *M;
	parse *p;

	p = parseCreateFileName(filename);

	if(!parseModelMicrostran(p,&M)){
		fprintf(stderr,"Failed to parse microstran model '%s'!\n", filename);
		exit(3);
	}

	parseDispose(p);

	fprintf(stderr,"Rendering frame with %d members...\n",M->num_membs);

	QWidget *mainwin = NULL;
	if(sceneoutfile){
		SoDB::init();
	}else{
		 mainwin = SoQt::init(argc, argv, argv[0]);
	}


	SoSeparator *root = new SoSeparator;

	root->addChild(axes());

	// render all the nodes in the model
	SbVec3f labeloffset(0.1,0.1,0.1);
	if(infotext){
		for(unsigned i=0; i<M->num_nodes; ++i){
			node_stmt *n;
			n = &(M->node[i]);
			stringstream ss;
			ss << n->id;
			SbVec3f p = labeloffset + vec3_to_coin(n->pos);
			//fprintf(stderr,"Label = %s\n",ss.str().c_str());
			//cerr << "pos = " << p << endl;
			root->addChild(text(p,ss.str().c_str(),RED));
		}
	}

	set<string> unknownsections;

	// render all the members in the model
	for(unsigned i=0; i<M->num_membs; ++i){
		memb_stmt *m;
		m = &(M->memb[i]);

		if(m->id < minmemb || m->id > maxmemb)continue;

		node_stmt *A = &(M->node[m->fromnode]);
		node_stmt *B = &(M->node[m->tonode]);
		prop_stmt *p = model_find_prop(M, m->prop);
	
		vec3 X = memb_get_orientation(M,m);

		SbVec3f vA = vec3_to_coin(A->pos);
		SbVec3f vB = vec3_to_coin(B->pos);
		SbVec3f vX = vec3_to_coin(X);
		vec3 dA, dB;

		vec3 ABn = vec3_norm(vec3_diff(B->pos,A->pos));

		moff_stmt *moff = NULL;
		if(memberoffsets){
			moff = model_find_member_offset(M, m->id);
			if(moff){
				if(moff->coordsys == MSTRANP_COORDS_LOCAL){
					ctrans_matrix c = ctrans_rotation_axes(ABn,X);

					if(A->id==5 && B->id==6){
						fprintf(stderr,"Coordinate transform:\n");
						ctrans_print(stderr,&c);
					}
					dA = ctrans_apply(c, moff->deltafrom);
					dB = ctrans_apply(c, moff->deltato);

					vec3 Y = vec3_norm(vec3_cross(ABn,X));
#if 0
					// work out the global coordinates another way...
					vec3 da1 = 
					 	vec3_add(
							vec3_scale(X, moff->deltafrom.x)
							,vec3_add(
								vec3_scale(Y, moff->deltafrom.y)
								,vec3_scale(ABn, moff->deltafrom.z)
							)
						);

					if(!vec3_equal_tol(da1,dA,1e-8)){
						VEC3_PR(ABn);
						VEC3_PR(X);
						VEC3_PR(Y);
						VEC3_PR(moff->deltafrom);
						VEC3_PR(da1);
						VEC3_PR(dA);
						CTRANS_PR(c);
						throw runtime_error("failed coordinate transformation");
					}
#endif
				}else{
					dA = moff->deltafrom;
					dB = moff->deltato;
				}
				vA += vec3_to_coin(dA);
				vB += vec3_to_coin(dB);
			}
		}
#if 0			
		if(model_get_member_offset_global(M, m->id, &moff)){
			vA += vec3_to_coin(moff.deltafrom);
			vB += vec3_to_coin(moff.deltato);
			cerr << "Member " << m->id << ": A offset = " << moff.deltafrom.x << "," << moff.deltafrom.y << "," << moff.deltafrom.z << endl;
		}
#endif

		//cerr << "Member " << m->id << " oriented to " << vX[0] << "," << vX[1] << "," << vX[2] << endl;

		stringstream ss;
		SbColor c;
		ss << m->id;
		if(p){
			const section *s = section_find(l, p->name);
			if(s){
				if(section_is_chs(s)){
					c = RED;
					double d = section_chs_outside_diameter(s) / 1000.; /* convert to metres */
					root->addChild(cylinder(vA,vB,d/2.,c));
				}else if(section_is_rod(s)){
					c = CYAN;
					root->addChild(cylinder(vA,vB,section_rod_diameter(s)/2.,c));
				}else if(section_is_isec(s)){
					c = GREEN;
					section_outline *o = section_isec_outline(s);
					root->addChild(prism(vA, vB, *o, c, vX));
					section_outline_destroy(o);
				}else if(section_is_shs(s)){
					c = YELLOW;
					section_outline *o = section_shs_outline(s);
					root->addChild(prism(vA, vB, *o, c, vX));
					section_outline_destroy(o);
				}else if(section_is_tophat(s)){
					c = ORANGE;
#if 0
					if(moff && (A->id<=40)){
						stringstream ss;
						ss << "dA global: (" << dA.x << "," << dA.y << "," << dA.z << ")" << endl;
						ss << "   local:  (" << moff->deltafrom.x << "," << moff->deltafrom.y << "," << moff->deltafrom.z << ")";

						root->addChild(arrow(vA - vec3_to_coin(vec3_norm(dA)), vA, CYAN, ss.str().c_str()));
						//root->addChild(arrow(vB - vec3_to_coin(vec3_norm(dB)), vB, PURPLE, "B"));
						root->addChild(arrow(vA - vec3_to_coin(dA), vA - vec3_to_coin(dA) + vec3_to_coin(X), RED, "X"));
						vec3 Y = vec3_norm(vec3_cross(ABn,X));
						root->addChild(arrow(vA - vec3_to_coin(dA), vA - vec3_to_coin(dA) + vec3_to_coin(Y), GREEN, "Y"));
						root->addChild(arrow(vA - vec3_to_coin(dA), vA - vec3_to_coin(dA) + vec3_to_coin(ABn), BLUE, "Z"));

						vec3 da1 = 
						 	vec3_add(
								vec3_scale(X, moff->deltafrom.x)
								,vec3_add(
									vec3_scale(Y, moff->deltafrom.y)
									,vec3_scale(ABn, moff->deltafrom.z)
								)
							);
						root->addChild(arrow(vA - vec3_to_coin(dA), vA - vec3_to_coin(dA) + vec3_to_coin(vec3_norm(da1)), PURPLE, "dA2"));
							
					}
#endif
					//cerr << "Rendering prism member" << endl;
					section_outline *o = section_tophat_outline(s);
					root->addChild(prism(vA, vB, *o, c, vX));
					section_outline_destroy(o);
				}
			}else{
				if(unknownsections.find(p->name)==unknownsections.end()){
					fprintf(stderr,"Warning: unknown section name '%s'\n",p->name);
					unknownsections.insert(p->name);
				}
				c = WHITE;
				root->addChild(cylinder(vA,vB,0.01/*radius*/,c));
			}
			ss << " (" << p->name << ")" << endl;
		}
		if(infotext){
			root->addChild(text(float(0.5)*(vA+vB)+labeloffset, ss.str().c_str(), c));
		}
	}

	root->ref();

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
		cerr << "Rendering scene to OpenGL directly..." << endl;

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



}
