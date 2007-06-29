#ifndef FRAME_H
#define FRAME_H

/**
	This is a 'frame' class that will hold a structure of
	nodes and edges. Initially we just want to be able to plot what
	it looks like, so it just needs to hold a group of point-like things
	plus a group of edge-like things, and track the relationships between
	them.
*/
class Frame{
public:
	Frame();
	Frame(const Frame &);

public:
	std::vector<Node> points;




#endif
