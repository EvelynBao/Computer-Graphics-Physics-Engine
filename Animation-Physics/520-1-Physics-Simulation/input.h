/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/


#ifndef _INPUT_H_
#define _INPUT_H_

//struct point userForce = { 0.0, 0.0, 0.0 };  // user force

/* Write a screenshot to the specified filename, in PPM format */
void saveScreenshot (int windowWidth, int windowHeight, char *filename);
void mouseButton(int button, int state, int x, int y);
void getMouseClickInfo(int* x, int* y, GLdouble* modelview, GLdouble* projection, GLint* view);

// mouse & keyboard control
void mouseMotionDrag(int x, int y);
void mouseMotion (int x, int y);
void mouseButton(int button, int state, int x, int y);
void keyboardFunc (unsigned char key, int x, int y);
struct point getMouseForce();
// read/write world files
void readWorld (char * fileName, struct world * jello);
void writeWorld (char * fileName, struct world * jello);

#endif

