/**@FILE
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

#include <sstream>
#include <set>
#include <string>
#include <stdexcept>
#include <iostream>
using namespace std;

#include "render.h"

#include <microstran/modelparser.h>
#include <microstran/sectionsparser.h>

const char *defaultsceneoutfile = "microstranmodel.iv";
const char *defaultlibfile = "microstran/properties.txt";

void usage(const char *progname){
	fprintf(stderr,"Usage: %s [-o [OUTFILE]] [-h] [-l LIBFILE] [-t] INFILE\n",progname);
	fprintf(stderr,"Load and parse a microstran file and render using Coin3D/OpenGL\n");
	fprintf(stderr,"  -h      Render using higher quality graphics (slower).\n");
	fprintf(stderr,"  INFILE  Microstran .arc file to render (eg 'model.arc')\n");
	fprintf(stderr,"  -l      Load a section library from LIBFILE.\n");
	fprintf(stderr,"  -t      Include text for node/member IDs and member sizes.\n");
	fprintf(stderr,"  OUTFILE Open Inventor file to output (defaults to '%s')\n\n",defaultsceneoutfile);
}

static SbVec3f vec3_to_coin(vec3 A){
	return SbVec3f(A.x,A.y,A.z);
}

int main(int argc, char **argv){

	const char *sceneoutfile = NULL;
	const char *libfile = defaultlibfile;
	bool highquality = false;
	bool infotext = false;

	// no options yet
	char c;
	while((c=getopt(argc,argv,"ho::l:t"))!=-1){
		switch(c){
			case 'h':
				highquality = 1;
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
	p1 = parseCreateFileName(libfile);
	section_library *l = NULL;
	l = section_library_create();
	parseSections(p1,l);
	fprintf(stderr,"done!\n");

	model *M;
	parse *p;

	p = parseCreateFileName(filename);

	if(!parseModelMicrostran(p,&M)){
		fprintf(stderr,"Failed to parse microstran model!\n");
		exit(3);
	}

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
		node_stmt *A = &(M->node[m->fromnode]);
		node_stmt *B = &(M->node[m->tonode]);
		prop_stmt *p = model_find_prop(M, m->prop);

		SbVec3f vA = vec3_to_coin(A->pos);
		SbVec3f vB = vec3_to_coin(B->pos);
		SbVec3f vX = vec3_to_coin(memb_get_orientation(M,m));

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
				}else if(section_is_isec(s)){
					c = GREEN;
					/* FIXME need to get the member orientation correct! */
					section_outline *o = section_isec_outline(s);
					cerr << "Member " << m->id << " oriented to " << vX[0] << "," << vX[1] << "," << vX[2] << endl;
					root->addChild(prism(vA, vB, *o, c, vX));
				}else if(section_is_shs(s)){
					c = YELLOW;
					/* FIXME need to get the member orientation correct! */
					section_outline *o = section_shs_outline(s);
					root->addChild(prism(vA, vB, *o, c, vX));
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
