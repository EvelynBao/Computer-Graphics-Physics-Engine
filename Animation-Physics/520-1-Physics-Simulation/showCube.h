/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/


#ifndef _SHOWCUBE_H_
#define _SHOWCUBE_H_

void showCube(struct world * jello);

void showBoundingBox();

void showInclinedPlane(struct world* jello);
struct point solvePlaneIntersection(double a, double b, double c, double d, double x, double y);

void showSkybox1(float size);
#endif
