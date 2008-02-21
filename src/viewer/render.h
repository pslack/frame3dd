#ifndef FRAME_RENDER_H
#define FRAME_RENDER_H

#include <microstran/sections.h>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/SbVec3f.h>

#include <vector>

extern const SbColor WHITE;
extern const SbColor RED;
extern const SbColor GREEN;
extern const SbColor YELLOW;
extern const SbColor BLUE;
extern const SbColor ORANGE;
extern const SbColor PURPLE;
extern const SbColor CYAN;

SoSeparator *text(const SbVec3f &left, const char *str, const SbColor &c=WHITE);
SoSeparator *sphere(const SbVec3f &C, const double &r=0.1, const SbColor &c=YELLOW);
SoSeparator *cylinder(const SbVec3f &A, const SbVec3f &B, const double &radius=0.1, const SbColor &c=YELLOW);
SoSeparator *cone(const SbVec3f &A, const SbVec3f &B, const double &r, const SbColor &c = YELLOW);
SoSeparator *axes(const double &size=1.0, double thickness=0.0, bool labelled=true);
SoSeparator *arrow(const SbVec3f &A, const SbVec3f &B, const SbColor &c=YELLOW, const char *label=NULL, double thickness=0.05);
SoSeparator *face(const SbVec3f &n, const std::vector<SbVec3f> &vertices, const SbColor &c=YELLOW);

SoSeparator *prism(const SbVec3f &A, const SbVec3f &B, const section_outline &outline, const SbColor &c, const SbVec3f &O);

#endif

