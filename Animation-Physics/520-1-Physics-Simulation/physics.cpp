/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include "input.h"
#include <algorithm>
#define NOMINMAX   
#include <windows.h>
/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */


void computeA(int p, int i, int j, int k, int neighbors[6][3], 
    struct world* jello, struct point p1, struct point v1, 
    double initialLength, double kElastic, double dElastic, struct point a[8][8][8]) {

    int ni = i + neighbors[p][0];
    int nj = j + neighbors[p][1];
    int nk = k + neighbors[p][2];

    // ignore those out of edge neighbours
    if (ni < 0 || ni >= 8 || nj < 0 || nj >= 8 || nk < 0 || nk >= 8)
        return;

    struct point p2 = jello->p[ni][nj][nk];
    struct point v2 = jello->v[ni][nj][nk];

    // position difference
    struct point pDiff;
    pDiff.x = p2.x - p1.x;
    pDiff.y = p2.y - p1.y;
    pDiff.z = p2.z - p1.z;

    // current length of string
    double L = sqrt(pDiff.x * pDiff.x + pDiff.y * pDiff.y + pDiff.z * pDiff.z);
    //if (i == 7 && j == 7 && k == 7) {
    //    printf("L at (%d,%d,%d) -> (%d,%d,%d) = %lf\n",
    //        i, j, k, ni, nj, nk, L);
    //}
    if (L == 0) return;
    //if (i == 7 && j == 7 && k == 7) {
    //    printf("L at (%d, %d, %d) -> (%d, %d, %d): %lf\n",
    //        i, j, k, ni, nj, nk, L);
    //}
    struct point unitDirection;
    unitDirection.x = pDiff.x / L;
    unitDirection.y = pDiff.y / L;
    unitDirection.z = pDiff.z / L;

    // F_s = -k * (L - L0) * unitDirection
    struct point F_s;
    double stretch = kElastic * (L - initialLength);
    //if (i == 7 && j == 7 && k == 7) {
    //    printf("stretch = kElastic * (L - initialLength): %lf\n",
    //        stretch);
    //}
    F_s.x = stretch * unitDirection.x;
    F_s.y = stretch * unitDirection.y;
    F_s.z = stretch * unitDirection.z;

    // velocity difference
    struct point vDiff;
    
    vDiff.x = v1.x - v2.x;
    vDiff.y = v1.y - v2.y;
    vDiff.z = v1.z - v2.z;
    
    // v dot L / L
    double vDotL = vDiff.x * unitDirection.x + vDiff.y * unitDirection.y + vDiff.z * unitDirection.z;

    //  F_d = -dElastic * v dot L / L * unitDirection
    struct point F_d;
    F_d.x = -dElastic * vDotL * unitDirection.x;
    F_d.y = -dElastic * vDotL * unitDirection.y;
    F_d.z = -dElastic * vDotL * unitDirection.z;

    // F = F_s + F_d
    struct point F;
    F.x = F_s.x + F_d.x;
    F.y = F_s.y + F_d.y;
    F.z = F_s.z + F_d.z;

    // a = F / m
    double mass = jello->mass;
    a[i][j][k].x += F.x / mass;
    a[i][j][k].y += F.y / mass;
    a[i][j][k].z += F.z / mass;
}


void computeStrcturalSprings(struct world* jello, struct point a[8][8][8]) {
    double initialLength = 1.0 / 7;
    double kElastic = jello->kElastic;
    double dElastic = jello->dElastic;

    int i, j, k;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 8; k++) {

                struct point p1 = jello->p[i][j][k]; // position
                struct point v1 = jello->v[i][j][k]; // velocity

                // 6 neighbours
                int neighbors[6][3] = {
                    {1, 0, 0}, {-1, 0, 0},
                    {0, 1, 0}, {0, -1, 0},
                    {0, 0, 1}, {0, 0, -1}
                };

                // for each neighbour of p1
                for (int p = 0; p < 6; p++) {
                    computeA(p, i, j, k, neighbors, jello, p1, v1, initialLength, kElastic, dElastic, a);
                }
            }
        }
    }
}

void computeShearSprings1(struct world* jello, struct point a[8][8][8]) {
    double initialLength = sqrt(2.0) / 7; 
    double kElastic = jello->kElastic;  
    double dElastic = jello->dElastic; 

    int i, j, k;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 8; k++) {

                struct point p1 = jello->p[i][j][k];
                struct point v1 = jello->v[i][j][k]; 

                // 12 neighbors
                int neighbors[12][3] = {
                    {1, 1, 0}, {-1, -1, 0}, {1, -1, 0}, {-1, 1, 0},
                    {0, 1, 1}, {0, -1, -1}, {0, 1, -1}, {0, -1, 1},
                    {1, 0, 1}, {-1, 0, -1}, {1, 0, -1}, {-1, 0, 1}
                };

                for (int n = 0; n < 12; n++) {
                    computeA(n, i, j, k, neighbors, jello, p1, v1, initialLength, kElastic, dElastic, a);
                }
            }
        }
    }
}

void computeShearSprings2(struct world* jello, struct point a[8][8][8]) {
    double initialLength = sqrt(3.0) / 7;
    double kElastic = jello->kElastic;
    double dElastic = jello->dElastic;

    int i, j, k;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 8; k++) {

                struct point p1 = jello->p[i][j][k];
                struct point v1 = jello->v[i][j][k]; 

                int neighbors[8][3] = {
                    {1, 1, 1}, {-1, -1, -1}, {1, -1, -1}, {-1, 1, 1},
                    {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}
                };

                for (int n = 0; n < 8; n++) {
                    computeA(n, i, j, k, neighbors, jello, p1, v1, initialLength, kElastic, dElastic, a);
                }
            }
        }
    }
}

void computeBendSprings(struct world* jello, struct point a[8][8][8]) {
    double initialLength = 2.0 / 7;    
    double kElastic = jello->kElastic; 
    double dElastic = jello->dElastic; 

    int i, j, k;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 8; k++) {

                struct point p1 = jello->p[i][j][k]; 
                struct point v1 = jello->v[i][j][k]; 

                // 6 neightbors
                int neighbors[6][3] = {
                    {2, 0, 0}, {-2, 0, 0},
                    {0, 2, 0}, {0, -2, 0},
                    {0, 0, 2}, {0, 0, -2}
                };

                for (int n = 0; n < 6; n++) {
                    computeA(n, i, j, k, neighbors, jello, p1, v1, initialLength, kElastic, dElastic, a);
                }
            }
        }
    }
}

// penalty acceleration
void computePA(int i, int j, int k, struct world* jello, struct point p1, struct point v1, struct point p2, 
    double initialLength, double kElastic, double dElastic, struct point normal, struct point a[8][8][8]) {

    struct point v2;
    v2.x = 0;
    v2.y = 0;
    v2.z = 0;

    // position difference
    struct point pDiff;
    pDiff.x = p2.x - p1.x;
    pDiff.y = p2.y - p1.y;
    pDiff.z = p2.z - p1.z;

    // current length of string
    double L = sqrt(pDiff.x * pDiff.x + pDiff.y * pDiff.y + pDiff.z * pDiff.z);

    if (L < 1e-6) return;  

    // unitDirection = direction / L
    struct point unitDirection;
    unitDirection.x = pDiff.x / L;
    unitDirection.y = pDiff.y / L;
    unitDirection.z = pDiff.z / L;

    unitDirection = normal;

    //  F_s = -k * (L - L0) * unitDirection
    struct point F_s;
    double stretch = kElastic * (L - initialLength);

    F_s.x = stretch * unitDirection.x;
    F_s.y = stretch * unitDirection.y;
    F_s.z = stretch * unitDirection.z;

    // velocity difference
    struct point vDiff;

    vDiff.x = v1.x - v2.x;
    vDiff.y = v1.y - v2.y;
    vDiff.z = v1.z - v2.z;

    // v dot L / L
    double vDotL = vDiff.x * unitDirection.x + vDiff.y * unitDirection.y + vDiff.z * unitDirection.z;

    //  F_d = -dElastic * v dot L / L * unitDirection
    struct point F_d;
    F_d.x = -dElastic * vDotL * unitDirection.x;
    F_d.y = -dElastic * vDotL * unitDirection.y;
    F_d.z = -dElastic * vDotL * unitDirection.z;

    // F = F_s + F_d
    struct point F;
    F.x = F_s.x + F_d.x;
    F.y = F_s.y + F_d.y;
    F.z = F_s.z + F_d.z;

    // a = F / m
    double mass = jello->mass;
    a[i][j][k].x += F.x / mass;
    a[i][j][k].y += F.y / mass;
    a[i][j][k].z += F.z / mass;
}

void checkCollision(struct world* jello, struct point a[8][8][8]) {

    double initialLength = 0;    // collision spring initial length
    double kElastic = jello->kElastic; 
    double dElastic = jello->dElastic;

    int i, j, k;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 8; k++) {

                struct point p1 = jello->p[i][j][k]; // current point position

                if (p1.x < -2 || p1.x>2 || p1.y < -2 || p1.y>2 || p1.z < -2 || p1.z>2) {

                    //struct point p1 = jello->p[i][j][k]; 
                    struct point v1 = jello->v[i][j][k]; 

                    struct point v2;
                    v2.x = 0;
                    v2.y = 0;
                    v2.z = 0;

                    if (p1.x < -2) {
                        struct point normal;
                        normal.x = 1;
                        normal.y = 0;
                        normal.z = 0;
                        struct point p2 = jello->p[i][j][k];
                        p2.x = -2;
                        computePA(i, j, k, jello, p1, v1, p2, initialLength, kElastic, dElastic, normal, a);
                    }
                    if (p1.x > 2) {
                        struct point normal;
                        normal.x = -1;
                        normal.y = 0;
                        normal.z = 0;
                        struct point p2 = jello->p[i][j][k];
                        p2.x = 2;
                        computePA(i, j, k, jello, p1, v1, p2, initialLength, kElastic, dElastic, normal, a);
                    }
                    if (p1.y < -2) {
                        struct point normal;
                        normal.x = 0;
                        normal.y = 1;
                        normal.z = 0;
                        struct point p2 = jello->p[i][j][k];
                        p2.y = -2;
                        computePA(i, j, k, jello, p1, v1, p2, initialLength, kElastic, dElastic, normal, a);
                    }
                    if (p1.y > 2) {
                        struct point normal;
                        normal.x = 0;
                        normal.y = -1;
                        normal.z = 0;
                        struct point p2 = jello->p[i][j][k];
                        p2.y = 2;
                        computePA(i, j, k, jello, p1, v1, p2, initialLength, kElastic, dElastic, normal, a);
                    }
                    if (p1.z < -2) {
                        struct point normal;
                        normal.x = 0;
                        normal.y = 0;
                        normal.z = 1;
                        struct point p2 = jello->p[i][j][k];
                        p2.z = -2;
                        computePA(i, j, k, jello, p1, v1, p2, initialLength, kElastic, dElastic, normal, a);
                    }
                    if (p1.z > 2) {
                        struct point normal;
                        normal.x = 0;
                        normal.y = 0;
                        normal.z = -1;
                        struct point p2 = jello->p[i][j][k];
                        p2.z = 2;
                        computePA(i, j, k, jello, p1, v1, p2, initialLength, kElastic, dElastic, normal, a);
                    }

                }
                if (jello->incPlanePresent == 1) {// inclined plane
                    double aa = jello->a, bb = jello->b, cc = jello->c, dd = jello->d;
                    double planeVal = aa * p1.x + bb * p1.y + cc * p1.z + dd;

                    if (planeVal < 0) {
                        struct point normal = { aa, bb, cc }; // normal of inclined plane
                        double normLength = sqrt(aa * aa + bb * bb + cc * cc);
                        normal.x /= normLength;
                        normal.y /= normLength;
                        normal.z /= normLength;

                        struct point p2 = p1;
                        double distance = fabs(planeVal) / normLength;
                        p2.x += normal.x * distance;
                        p2.y += normal.y * distance;
                        p2.z += normal.z * distance;

                        computePA(i, j, k, jello, p1, jello->v[i][j][k], p2, initialLength, kElastic, dElastic, normal, a);
                    }              
                }
            }
        }
    }
}

void computeForceField(struct world* jello, struct point a[8][8][8]) {
    int res = jello->resolution;
    double mass = jello->mass;

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                // get position
                struct point p = jello->p[i][j][k];

                // normalize to [0, resolution-1]
                double u = (p.x + 2) / 4.0 * (res - 1);
                double v = (p.y + 2) / 4.0 * (res - 1);
                double w = (p.z + 2) / 4.0 * (res - 1);

                // 8 nearest force points
                int i0 = floor(u), i1 = i0 + 1;
                int j0 = floor(v), j1 = j0 + 1;
                int k0 = floor(w), k1 = k0 + 1;

                // limit range
                i0 = (i0 < 0) ? 0 : ((i0 > res - 1) ? (res - 1) : i0);
                i1 = (i1 < 0) ? 0 : ((i1 > res - 1) ? (res - 1) : i1);
                j0 = (j0 < 0) ? 0 : ((j0 > res - 1) ? (res - 1) : j0);
                j1 = (j1 < 0) ? 0 : ((j1 > res - 1) ? (res - 1) : j1);
                k0 = (k0 < 0) ? 0 : ((k0 > res - 1) ? (res - 1) : k0);
                k1 = (k1 < 0) ? 0 : ((k1 > res - 1) ? (res - 1) : k1);

                // weight
                double du = u - i0;
                double dv = v - j0;
                double dw = w - k0;

                // force point value
                struct point f000 = jello->forceField[i0 * res * res + j0 * res + k0];
                struct point f100 = jello->forceField[i1 * res * res + j0 * res + k0];
                struct point f010 = jello->forceField[i0 * res * res + j1 * res + k0];
                struct point f110 = jello->forceField[i1 * res * res + j1 * res + k0];
                struct point f001 = jello->forceField[i0 * res * res + j0 * res + k1];
                struct point f101 = jello->forceField[i1 * res * res + j0 * res + k1];
                struct point f011 = jello->forceField[i0 * res * res + j1 * res + k1];
                struct point f111 = jello->forceField[i1 * res * res + j1 * res + k1];

                // calculate Fx, Fy, Fz
                struct point F;
                F.x = f000.x * (1 - du) * (1 - dv) * (1 - dw) +
                      f100.x * du * (1 - dv) * (1 - dw) +
                      f010.x * (1 - du) * dv * (1 - dw) +
                      f110.x * du * dv * (1 - dw) +
                      f001.x * (1 - du) * (1 - dv) * dw +
                      f101.x * du * (1 - dv) * dw +
                      f011.x * (1 - du) * dv * dw +
                      f111.x * du * dv * dw;

                F.y = f000.y * (1 - du) * (1 - dv) * (1 - dw) +
                      f100.y * du * (1 - dv) * (1 - dw) +
                      f010.y * (1 - du) * dv * (1 - dw) +
                      f110.y * du * dv * (1 - dw) +
                      f001.y * (1 - du) * (1 - dv) * dw +
                      f101.y * du * (1 - dv) * dw +
                      f011.y * (1 - du) * dv * dw +
                      f111.y * du * dv * dw;

                F.z = f000.z * (1 - du) * (1 - dv) * (1 - dw) +
                      f100.z * du * (1 - dv) * (1 - dw) +
                      f010.z * (1 - du) * dv * (1 - dw) +
                      f110.z * du * dv * (1 - dw) +
                      f001.z * (1 - du) * (1 - dv) * dw +
                      f101.z * du * (1 - dv) * dw +
                      f011.z * (1 - du) * dv * dw +
                      f111.z * du * dv * dw;

                // a = F / m
                a[i][j][k].x += F.x / mass;
                a[i][j][k].y += F.y / mass;
                a[i][j][k].z += F.z / mass;
            }
        }
    }
    

}

void computeMouseForce(struct world* jello, struct point a[8][8][8]) {
    // ignore all non left mouse
    if (!g_iLeftMouseButton) {
        return;
    }
    struct point mouseForce = getMouseForce();
    if (mouseForce.x == 0 && mouseForce.y == 0 && mouseForce.z == 0) {
        return;
    }
    double mass = jello->mass;

    struct point worldForce;
    worldForce.x = 0;
    worldForce.y = mouseForce.x;
    worldForce.z = mouseForce.y;

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                a[i][j][k].x += worldForce.x / mass;
                a[i][j][k].y += worldForce.y / mass;
                a[i][j][k].z += worldForce.z / mass;
            }
        }
    }

    // other drag method

    /*
    //int mouseX, mouseY;
    //GLdouble modelview[16], projection[16];
    //GLint viewport[4];

    //getMouseClickInfo(&mouseX, &mouseY, modelview, projection, viewport);

    // read depth
    //GLfloat depth;
    //glReadPixels(mouseX, mouseY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

    // calculate world position
    //GLdouble worldX, worldY, worldZ;
    //gluUnProject(mouseX, mouseY, depth, modelview, projection, viewport, &worldX, &worldY, &worldZ);

    //printf("[computeMouseForce] Clicked 2D: (%d, %d), Depth: %f → 3D: (%.2lf, %.2lf, %.2lf)\n",
    //    mouseX, mouseY, depth, worldX, worldY, worldZ);

    // find the nearest point on jello cube
    //struct point nearestPoint;
    //double minDist = 99999.0;
    */
    
    /*
    // apply mouse force to the central 8 points
    //int centerIndices[8][3] = {
    //    {3, 3, 3}, {3, 3, 4}, {3, 4, 3}, {3, 4, 4},
    //    {4, 3, 3}, {4, 3, 4}, {4, 4, 3}, {4, 4, 4}
    //};

    
    int centerIndices[8][3] = {
    {0, 0, 0}, {0, 7, 0}, {7, 0, 0}, {7, 7, 0},
    {0, 0, 7}, {0, 7, 7}, {7, 0, 7}, {7, 7, 7}
    };

    for (int n = 0; n < 8; n++) {
        int i = centerIndices[n][0];
        int j = centerIndices[n][1];
        int k = centerIndices[n][2];

        a[i][j][k].x += worldForce.x / mass;
        a[i][j][k].y += worldForce.y / mass;
        a[i][j][k].z += worldForce.z / mass;

        printf("Force applied to (%d, %d, %d): (%.2f, %.2f, %.2f)\n", i, j, k,
            worldForce.x / mass, worldForce.y / mass, worldForce.z / mass);
        printf("→ Acceleration at (%d, %d, %d): a = (%.4f, %.4f, %.4f)\n",
            i, j, k, a[i][j][k].x, a[i][j][k].y, a[i][j][k].z);
        printf("[DEBUG] v[%d][%d][%d] = (%.4lf, %.4lf, %.4lf)\n", i, j, k,
            jello->v[i][j][k].x, jello->v[i][j][k].y, jello->v[i][j][k].z);
        printf("[DEBUG] p[%d][%d][%d] = (%.4lf, %.4lf, %.4lf)\n", i, j, k,
            jello->p[i][j][k].x, jello->p[i][j][k].y, jello->p[i][j][k].z);
    }
    
    //printf("User Force: x = %lf, y = %lf, z = %lf\n", mouseForce.x, mouseForce.y, mouseForce.z);
*/
 
    /* no 2d -> 3d
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                a[i][j][k].x += mouseForce.x / mass;
                a[i][j][k].y += mouseForce.y / mass;
                a[i][j][k].z += mouseForce.z / mass;
            }
        }
    }
    */

    /* apply force on the surface of the cube(nearest point)
    int ni = 0, nj = 0, nk = 0;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                struct point p = jello->p[i][j][k];

                double dx = p.x - worldX;
                double dy = p.y - worldY;
                double dz = p.z - worldZ;
                double dist = sqrt(dx * dx + dy * dy + dz * dz);

                if (dist < minDist) {
                    minDist = dist;
                    nearestPoint = p;
                    ni = i; nj = j; nk = k;
                }
            }
        }
    }
    //printf("[computeMouseForce] Nearest point in jello: (%d, %d, %d) → (%.2lf, %.2lf, %.2lf)\n",
    //    ni, nj, nk, nearestPoint.x, nearestPoint.y, nearestPoint.z);

    double influenceRadius = 0.4;
    struct point worldForce;
    worldForce.x = mouseForce.x * modelview[0] + mouseForce.y * modelview[4];
    worldForce.y = mouseForce.x * modelview[1] + mouseForce.y * modelview[5];
    worldForce.z = mouseForce.x * modelview[2] + mouseForce.y * modelview[6];
    printf("World Force: x=%.2f, y=%.2f, z=%.2f\n", worldForce.x, worldForce.y, worldForce.z);

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                struct point p = jello->p[i][j][k];

                double dx = p.x - nearestPoint.x;
                double dy = p.y - nearestPoint.y;
                double dz = p.z - nearestPoint.z;
                double distance = sqrt(dx * dx + dy * dy + dz * dz);

                double weight = exp(-distance / influenceRadius);

                a[i][j][k].x += weight * worldForce.x / mass;
                a[i][j][k].y += weight * worldForce.y / mass;
                a[i][j][k].z += weight * worldForce.z / mass;
            }
        }
    }
    */

}

void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
    int i, j, k;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 8; k++) {
                a[i][j][k].x = 0; 
                a[i][j][k].y = 0;
                a[i][j][k].z = 0;
            }
        }
    }

    computeStrcturalSprings(jello, a);
    computeShearSprings1(jello, a);
    computeShearSprings2(jello, a);
    computeBendSprings(jello, a);

    checkCollision(jello, a);

    computeForceField(jello, a);

    computeMouseForce(jello, a);

    //printf("Acceleration at (4,4,4): %lf, %lf, %lf\n",
    //    a[4][4][4].x, a[4][4][4].y, a[4][4][4].z);

}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
