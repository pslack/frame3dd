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
*//* @FILE
	Some rendering routines for use with the Open Inventor API that allow
	high-level objects to be added to the scenegraph easily.
*/
#include "render.h"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoCone.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoLOD.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoNormalBinding.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoShapeHints.h>

extern "C"{
#include <microstran/array.h>
#include <microstran/vec3.h>
};

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
using namespace std;

const SbColor WHITE(1,1,1);
const SbColor RED(1,0,0);
const SbColor GREEN(0,1,0);
const SbColor YELLOW(1,1,0);
const SbColor BLUE(0.5,0.5,1);
const SbColor ORANGE(1,0.6,0);
const SbColor PURPLE(160./255,32./255,240./255);
const SbColor MAGENTA(1,0,1);
const SbColor CYAN(0,1,1);
const SbColor LIME(0.71,1,0);

static float angle(const SbVec3f &a, const SbVec3f &b){
	SbVec3f c(a.cross(b));
	return atan2(c.length(),a.dot(b));
}

SbVec3f vec3_to_coin(vec3 A){
	return SbVec3f(A.x,A.y,A.z);
}

static vec3 vec3_from_coin(SbVec3f A){
	return vec3_create(A[0],A[1],A[2]);
}

SoSeparator *text(const SbVec3f &left, const char *str, const SbColor &c){
	SoSeparator *s = new SoSeparator;

	SoBaseColor *col = new SoBaseColor;
	col->rgb = c;
	s->addChild(col);

	SoTranslation *tr = new SoTranslation;
	tr->translation = left;
	s->addChild(tr);

	SoText2 *txt = new SoText2;
	txt->string.setValue(str);
	s->addChild(txt);

	return s;
}

SoSeparator *sphere(const SbVec3f &C, const double &r, const SbColor &c){
	SoSeparator *s = new SoSeparator;

	SoBaseColor *col = new SoBaseColor;
	col->rgb = c;
	s->addChild(col);

	SoTranslation *tr = new SoTranslation;
	tr->translation = C;
	s->addChild(tr);

	SoSphere *sph = new SoSphere;
	sph->radius = r;
	s->addChild(sph);

	return s;
}

SoSeparator *cylinder(const SbVec3f &A, const SbVec3f &B, const double &radius, const SbColor &c){

	SoSeparator *s = new SoSeparator;

	SbVec3f AB(B - A);
	AB.normalize();
	SbVec3f d(AB);
	SbVec3f m(float(0.5) * (A + B));

	SbVec3f z(0,1,0);
	SbVec3f a = z.cross(d);
	double theta = angle(z,d);
	if(a.length() < 1e-6){
		a = SbVec3f(0,0,1);
		theta = 0;
	}

	SoBaseColor *col = new SoBaseColor;
	col->rgb = c;
	s->addChild(col);

	//cerr << "translation = " << m << endl;
	SoTranslation *tr = new SoTranslation;
	tr->translation = m;
	s->addChild(tr);

	//cerr << "rotation = " << theta * 180./PI << " deg around " << a << endl;
	SoRotation *rot = new SoRotation;
	rot->rotation = SbRotation(a,theta);
	s->addChild(rot);

	//SoComplexity *cplx = new SoComplexity;
	//cplx->type = SoComplexity::SCREEN_SPACE;
	//s->addChild(cplx);

	SoCylinder *cyl = new SoCylinder;
	cyl->height = (A-B).length();
	cyl->radius = radius;
	s->addChild(cyl);

	return s;
}

/**
	Cone with its axis on AB, with its base at A and its apex at B. The base
	radius is r.
*/
SoSeparator *cone(const SbVec3f &A, const SbVec3f &B, const double &r, const SbColor &c){

	SoSeparator *s = new SoSeparator;

	SbVec3f AB(B - A);
	AB.normalize();
	SbVec3f d(AB);
	SbVec3f m(float(0.5) * (A + B));

	SbVec3f z(0,1,0);
	SbVec3f a = z.cross(d);
	double theta = angle(z,d);
	if(a.length() < 1e-6){
		a = SbVec3f(0,0,1);
		theta = 0;
	}

	SoBaseColor *col = new SoBaseColor;
	col->rgb = c;
	s->addChild(col);

	//cerr << "translation = " << m << endl;
	SoTranslation *tr = new SoTranslation;
	tr->translation = m;
	s->addChild(tr);

	//cerr << "rotation = " << theta * 180./PI << " deg around " << a << endl;
	SoRotation *rot = new SoRotation;
	rot->rotation = SbRotation(a,theta);
	s->addChild(rot);

	SoCone *cyl = new SoCone;
	cyl->height = (B-A).length();
	cyl->bottomRadius = r;
	s->addChild(cyl);

	return s;
}

SoSeparator *axes(const double &size, double thickness,bool labelled){

	if(thickness==0){
		thickness = size/20.;
	}

 	SoSeparator *s = new SoSeparator;
	s->addChild(cylinder(SbVec3f(0,0,0), SbVec3f(size,0,0), thickness, RED));
	s->addChild(cone(SbVec3f(size,0,0), SbVec3f(size+3*thickness,0,0), 2*thickness, RED));

	s->addChild(cylinder(SbVec3f(0,0,0), SbVec3f(0,size,0), thickness, GREEN));
	s->addChild(cone(SbVec3f(0,size,0), SbVec3f(0,size+3*thickness,0), 2*thickness, GREEN));

	s->addChild(cylinder(SbVec3f(0,0,0), SbVec3f(0,0,size), thickness, BLUE));
	s->addChild(cone(SbVec3f(0,0,size), SbVec3f(0,0,size+3*thickness), 2*thickness, BLUE));

	if(labelled)s->addChild(text(SbVec3f(0,0,size+5*thickness),"Z"));
	if(labelled)s->addChild(text(SbVec3f(0,size+5*thickness,0),"Y"));
	if(labelled)s->addChild(text(SbVec3f(size+5*thickness,0,0),"X"));

	return s;
}

SoSeparator *arrow(const SbVec3f &A, const SbVec3f &B, const SbColor &c, const char *label, double thickness){
	SoSeparator *s = new SoSeparator;
	SbVec3f Bs = B + float(3. * thickness) * (A-B);
	s->addChild(cylinder(A,Bs,0.5*thickness, c));
	s->addChild(cone(Bs,B,thickness, c));
	if(label!=NULL){
		s->addChild(text(float(0.5)*(A+Bs)+SbVec3f(thickness,thickness,thickness), label, c));
	}
	return s;
}

/**
	Return a 'face set' containing just a single face.
	@param n direction of the surface normal for the face
	@param vertices the vertices lying on the face
	@param c color in which to render it.
*/
SoSeparator *face(const SbVec3f &n, const vector<SbVec3f> &vertices, const SbColor &c){
	SoSeparator *s = new SoSeparator;

	SoNormal *n1 = new SoNormal;
	n1->vector.setValues(0,1, &n);
	s->addChild(n1);

	SoNormalBinding *nb = new SoNormalBinding;
	nb->value = SoNormalBinding::PER_FACE;
	s->addChild(nb);

	SoBaseColor *col = new SoBaseColor;
	col->rgb = c;
	s->addChild(col);

	SoCoordinate3 *cor = new SoCoordinate3;
	cor->point.setValues(0,vertices.size(),&(vertices[0]));
	s->addChild(cor);

	SoFaceSet *fs = new SoFaceSet;
	const int32_t nv = vertices.size();
	fs->numVertices.setValues(0,1,&nv);
	s->addChild(fs);

	return s;
}

/**
	Render an arbitrary right-prism length

	@param A first endpoint
	@param B second endpoint
	@param o section outline, tells how to draw the section.
	@param O orientation vector for the +X axis of the local coordinate system (see memb_get_orientation for the way that this is calculated)

	@TODO FIXME the orientation code should be simplified using ctrans_rotation_axes
*/
SoSeparator *prism(const SbVec3f &A, const SbVec3f &B, const section_outline_struct &o, const SbColor &c, const SbVec3f &O){
	SoSeparator *s = new SoSeparator;

	SoBaseColor *col = new SoBaseColor;
	col->rgb = c;
	s->addChild(col);

	SbVec3f minax = O;
	minax.normalize();
	SbVec3f majax = (B-A).cross(O);
	majax.normalize();
	//s->addChild(arrow(A,A+majax,GREEN,"MAJOR"));
	//s->addChild(arrow(A,A+minax,RED,"MINOR"));

	// member half-length
	double r = 0.5*(B-A).length();
	unsigned n = ARRAY_NUM(o.trace);
	unsigned nv = ARRAY_NUM(o.point);

	// direction vector
	SbVec3f d(B - A);
	d.normalize();

	// midpoint of the member
	SbVec3f m(float(0.5) * (A + B));

	// as-placed orientation of cylinder
	SbVec3f y(0,0,1);

	// axis and angle of rotation to get the right position
	SbVec3f a = y.cross(d);
	double theta_cyl = angle(y,d);
	// fix case where cross product is 0
	if(a.length() < 1e-6){
		a = SbVec3f(1,0,0);
		theta_cyl = 0;
	}

	//s->addChild(arrow(SbVec3f(0,0,0),a,CYAN,"rotation"));

	vec3 oo = vec3_norm(vec3_from_coin(majax));
#if 0
	vec3 vz = vec3_create(0,0,1);
	vec3 vy = vec3_create(0,1,0);
	vec3 vx = vec3_create(1,0,0);
	assert(vec3_mod(vec3_diff(vec3_rotate(vy,vx,M_PI/2.),vz))< 1e-8);
	assert(vec3_mod(vec3_diff(vec3_rotate(vy,vx,-M_PI/2.),vec3_scale(vz,-1)))< 1e-8);
#endif

	// FIXME still need to fix the twist-angle stuff to get the member oriented correctly
	SbVec3f Ad_unrot = vec3_to_coin(vec3_rotate(oo,vec3_from_coin(a),-theta_cyl));

	//s->addChild(arrow(A,A+Ad_unrot,YELLOW,"unrotated major axis"));

	vec3 om = vec3_norm(vec3_from_coin(minax));
	SbVec3f om_unrot = vec3_to_coin(vec3_rotate(om,vec3_from_coin(a),-theta_cyl));
	//s->addChild(arrow(A,A+om_unrot,YELLOW,"unrotated minor axis"));

# if 0
	cerr << "dot = "<< Ad_unrot.dot(SbVec3f(1,0,0)) << endl;
	if(fabs(Ad_unrot.dot(SbVec3f(1,0,0)))>1e-8){
		cerr << "Unrotated = " << Ad_unrot[0] << ", " << Ad_unrot[1] << ", " << Ad_unrot[2] << endl;
		throw runtime_error("Orientation O is wrong for prism");
	}
# endif

	SbVec3f zeroA(0,1,0);
	double theta_A = angle(zeroA, Ad_unrot);
	SbVec3f dir_A = zeroA.cross(Ad_unrot);

	//s->addChild(arrow(A,A+zeroA,PURPLE,"zeroA"));

	//cerr << "theta_A = " << theta_A * 180/M_PI << endl;

	//SoLevelOfDetail here...
	SoLOD *lod = new SoLOD;
	float rvals[] = {60};
	lod->center = float(0.5)*(A+B);
	lod->range.setValues(0,1,rvals);

	// DETAILED level of detail
	SoSeparator *detailed = new SoSeparator;

#if ORIENTATION_DEBUG
	detailed->addChild(arrow(A, A+minax, RED, "X"));
	detailed->addChild(arrow(A, A+majax, GREEN, "Y"));
	detailed->addChild(arrow(A, A+d, BLUE, "Z"));
#endif

#if 1
	SoShapeHints *sha = new SoShapeHints;
	sha->vertexOrdering = SoShapeHints::UNKNOWN_ORDERING;
	sha->shapeType = SoShapeHints::UNKNOWN_SHAPE_TYPE;
    sha->faceType = SoShapeHints::UNKNOWN_FACE_TYPE;
	detailed->addChild(sha);
#endif

	// rotation to orient the prepared member correctly in space
	SoTransform *tr = new SoTransform;
	tr->translation.setValue(m);
	tr->rotation = SbRotation(a,theta_cyl);
	detailed->addChild(tr);

	// rotation to orient the prepared member correctly around its axis
	SoTransform *axi = new SoTransform;
	axi->rotation = SbRotation(dir_A,theta_A);
	detailed->addChild(axi);

	SbVec3f *v3;
	v3 = new SbVec3f[2*nv];
	for(unsigned i=0; i<nv; ++i){
		vec2 *v2 = (vec2 *)array_get((array *)&(o.point), i);
		// create the vertices at the +z end...
		v3[i] = SbVec3f(v2->x, v2->y, +r);
		// create the vertices at the -z end...
		v3[nv + i] = SbVec3f(v2->x,v2->y, -r);
	}

#if 0
	SoMaterialBinding *mbi = new SoMaterialBinding;
	mbi->value = SoMaterialBinding::PER_FACE;
	detailed->addChild(mbi);
#endif

	SoCoordinate3 *coo = new SoCoordinate3;
	coo->point.setValues(0, 2 * nv, v3);
	detailed->addChild(coo);

	unsigned sides = 0;
	for(unsigned i=0; i<n-1; ++i){
		const int &v1 = *((int *)array_get((array *)&(o.trace), i));
		if(v1==-1)continue; /* skip this one, if last one was a 'pen up' */
		const int &v2 = *((int *)array_get((array *)&(o.trace), i+1));
		if(v2==-1)continue; /* skip this one if we're about to pen-up*/
		sides++;
	}

	int *i3; /* indices for the 3D face set */
	i3 = new int[5*sides + 2 * (sides + 1)];

	// add the sides of the prism
	unsigned ni=0;
	for(unsigned i=0; i<n-1; ++i){
		int v1 = *(const int *)array_get((array *)&(o.trace), i);
		int v2 = *(const int *)array_get((array *)&(o.trace), i+1);
		if(v2==-1 || v1==-1)continue;
		//cerr << "Side from v" << v1 << " to v" << v2 << endl;
		i3[ni++] = v1;
		i3[ni++] = v2;
		i3[ni++] = nv + v2;
		i3[ni++] = nv + v1;
		i3[ni++] = SO_END_FACE_INDEX;
	}

#if 1
	SoShapeHints *sha1 = new SoShapeHints;
	sha1->vertexOrdering = SoShapeHints::CLOCKWISE;
	sha1->shapeType = SoShapeHints::UNKNOWN_SHAPE_TYPE;
    sha1->faceType = SoShapeHints::UNKNOWN_FACE_TYPE;
	detailed->addChild(sha1);
#endif

	// add the ends of the prism
#if 1
	for(unsigned i=0; i<n-1; ++i){
		const int &v1 = *((int *)array_get((array *)&(o.trace), i));
		if(v1==-1)continue; /* skip this one, if last one was a 'pen up' */
		const int &v2 = *((int *)array_get((array *)&(o.trace), i+1));
		if(v2==-1)continue; /* skip this one if we're about to pen-up*/
		i3[ni++] = v1;
		//cerr << "v" << v1 << "=" << coo->point[v1][0] << "," << coo->point[v1][1] << endl;
		//cerr << " to v" << v2 << "=" << coo->point[v2][0] << "," << coo->point[v2][1] << endl;
	}
	//cerr << endl;

	for(unsigned i=0; i<n-1; ++i){
		const int &v1 = *((int *)array_get((array *)&(o.trace), i));
		if(v1==-1)continue; /* skip this one, if last one was a 'pen up' */
		const int &v2 = *((int *)array_get((array *)&(o.trace), i+1));
		if(v2==-1)continue; /* skip this one if we're about to pen-up*/
		i3[ni++] = nv + v1;
		//cerr << "v" << v1 << "=" << coo->point[v1][0] << "," << coo->point[v1][1] << endl;
		//cerr << " to v" << v2 << "=" << coo->point[v2][0] << "," << coo->point[v2][1] << endl;
	}
	i3[ni++] = SO_END_FACE_INDEX;

#endif

	SoIndexedFaceSet *fsi = new SoIndexedFaceSet;
	fsi->coordIndex.setValues(0, ni, i3);
	detailed->addChild(fsi);
	lod->addChild(detailed);

	delete[] i3;
	delete[] v3;

	// rough model - just a cylinder?
	SoSeparator *rough = new SoSeparator;
	double d_rough = section_outline_approx_diameter(&o);
	rough->addChild(cylinder(A,B,d_rough/2.,c));
	lod->addChild(rough);
	s->addChild(lod);

	return s;
}


#ifdef RENDER_DEBUG
int
main(int argc, char ** argv){

    // Initializes SoQt library (and implicitly also the Coin and Qt
    // libraries). Returns a top-level / shell Qt window to use.
    QWidget *mainwin = SoQt::init(argc, argv, argv[0]);

    // Make a simple scene graph by using the Coin library
    SoSeparator *root = new SoSeparator;
    root->ref();

	root->addChild(axes(3,0.01));

	SoComplexity *cplx = new SoComplexity;
	cplx->type = SoComplexity::SCREEN_SPACE;
	root->addChild(cplx);

	{
#if 0
		Plane A(Vector(1,1,0),Vector(-2,1,0));
		Plane B(Vector(1,1,4),Vector(1,0,0));
#endif

#if 0
		srand(clock());

		Vector A(1,1,1);
		Vector dA(0,1,0);
		double rA = (double(rand())/RAND_MAX - 1./2) * 360;
		dA.rotate(Vector(-1,0,0),rA);
		Plane PA(A,normalise(dA));

		double rB = (double(rand())/RAND_MAX - 1./2) * 360;
		Vector B(3.5,1,1);
		Vector dB(0,1,0);
		dB.rotate(Vector(1,0,0),rB);
		Plane PB(B,normalise(dB));

		Vector ax(1,0,0);
		A.rotate(ax,45);
		B.rotate(ax,45);

		cerr << "A = " << PA;
		cerr << "B = " << PB;
		root->addChild(squashedpipe(PA,PB,0.2,YELLOW,0.004, 0.7));

		root->addChild(sphere(0.5*(PA.p + PB.p),0.15));
		root->addChild(sphere(PA.p,0.05,RED));
		root->addChild(sphere(PB.p,0.05,GREEN));

		root->addChild(arrow(PA.p,PA.p + PA.d,RED,"A dir"));
		root->addChild(arrow(PB.p,PB.p + PB.d,GREEN,"B dir"));
#endif
	}

#if 0
	{
		Plane B(Vector(0,0,0),Vector(0,0,1));
		Plane A(Vector(5,5,0),Vector(1./sqrt(2),-1./sqrt(2),1./sqrt(2)));
		root->addChild(squashedpipe(A,B,0.2,WHITE));
	}
#endif

    // Use one of the convenient SoQt viewer classes.
    SoQtExaminerViewer *eviewer = new SoQtExaminerViewer(mainwin);
    eviewer->setSceneGraph(root);
    eviewer->show();

    // Pop up the main window.
    SoQt::show(mainwin);
    // Loop until exit.
    SoQt::mainLoop();

    // Clean up resources.
    delete eviewer;
    root->unref();

    return 0;
}
#endif

