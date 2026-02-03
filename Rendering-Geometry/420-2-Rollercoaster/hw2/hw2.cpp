/*
  CSCI 420 Computer Graphics, USC
  Assignment 2: Roller Coaster
  C/C++ starter code

  Student username: ruyibao
*/

#include "openGLHeader.h"
#include "glutHeader.h"
#include "openGLMatrix.h"
#include "imageIO.h"
#include "pipelineProgram.h"
#include "vbo.h"
#include "vao.h"

#include <iostream>
#include <cstring>
#include <windows.h> // for sleep function
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openGLHeader.h"
#include "imageIO.h"
#include <vector>
#include <thread>
#include <future>
#include <vector>
#include <chrono> 
//#include <glm/vec3.hpp>
//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp> // For glm::translate, rotate, scale, perspective
//#include <glm/gtc/type_ptr.hpp> // For glm::value_ptr

#if defined(WIN32) || defined(_WIN32)
  #ifdef _DEBUG
    #pragma comment(lib, "glew32d.lib")
  #else
    #pragma comment(lib, "glew32.lib")
  #endif
#endif

#if defined(WIN32) || defined(_WIN32)
  char shaderBasePath[1024] = SHADER_BASE_PATH;
#else
  char shaderBasePath[1024] = "../openGLHelper";
#endif

using namespace std;

int mousePos[2]; // x,y screen coordinates of the current mouse position

int leftMouseButton = 0; // 1 if pressed, 0 if not 
int middleMouseButton = 0; // 1 if pressed, 0 if not
int rightMouseButton = 0; // 1 if pressed, 0 if not

bool takeScreenshot = false;

// Represents one spline control point.
struct Point
{
    double x, y, z;
};

// Contains the control points of the spline.
struct Spline
{
    int numControlPoints;
    Point* points;
} spline;


typedef enum { ROTATE, TRANSLATE, SCALE } CONTROL_STATE;
CONTROL_STATE controlState = ROTATE;

//typedef enum { ONE_POINT, TWO_LINE, THREE_TRIANGLE, FOUR_SMOOTH } RENDERING_MODE;
//RENDERING_MODE renderingMode;

typedef enum { NOTSMOOTH, SMOOTH } SHADER_MODE;
SHADER_MODE shaderMode;

// Transformations of the terrain.
float terrainRotate[3] = { 0.0f, 0.0f, 0.0f }; 
// terrainRotate[0] gives the rotation around x-axis (in degrees)
// terrainRotate[1] gives the rotation around y-axis (in degrees)
// terrainRotate[2] gives the rotation around z-axis (in degrees)
float terrainTranslate[3] = { 0.0f, 0.0f, 0.0f };
float terrainScale[3] = { 1.0f, 1.0f, 1.0f };

float scale = 1.0f;
float exponent = 1.0f;

// Width and height of the OpenGL window, in pixels.
int windowWidth = 1280;
int windowHeight = 720;
char windowTitle[512] = "CSCI 420 Homework 2";

// Stores the image loaded from disk.
ImageIO * heightmapImage;
int heightmapImageHeight = 0;
int heightmapImageWidth = 0;

// Stores the image loaded from disk.
ImageIO* groundImage;
int groundImageHeight = 0;
int groundImageWidth = 0;

// Stores the image loaded from disk.
ImageIO* skyImage;
int skyImageHeight = 0;
int skyImageWidth = 0;

// Number of vertices in the single triangle (starter code).
int numVertices;
int numVerticesGround;
int numVerticesSky;

// number of vertices in one tube
int numTubeVertices;
int numTube;


// CSCI 420 helper classes.
OpenGLMatrix matrix;
PipelineProgram * pipelineProgram = nullptr;
PipelineProgram* pipelineProgramGround = nullptr;
PipelineProgram* pipelineProgramSky = nullptr;


// rail
VBO * vboVertices1 = nullptr;
VBO * vboNormals1 = nullptr;
VAO * vao1 = nullptr;
//ground
VBO* vboVerticesGround = nullptr;
VBO* vboColorsGround = nullptr;
VAO* vaoGround = nullptr;
//sky
VBO* vboVerticesSky = nullptr;
VBO* vboColorsSky = nullptr;
VAO* vaoSky = nullptr;

double cameraSpeed = 0.02;
double moveTime = 0.06; // 16.666 FPS
double eyez = 15.0;
double p = 0.0;
double gravity = 9.8;                                                             
double hmax = 0.0;

float unit = 0.1f;


// Track Vertex positions and normals.
float* positions;
float* colortonormals;

// Ground Vertex positions and normals.
float* positionsGround;
float* colorsGround;

// Sky Ground Vertex positions and normals.
float* positionsSky;
float* colorsSky;

vector<Point> splinePoints;
vector<Point> splineTangents;
vector<Point> splineNormalizedTangents;
vector<Point> splineNormals;
vector<Point> splineBinormals;

GLuint textureHandle;

// Write a screenshot to the specified filename.
void saveScreenshot(const char * filename)
{
    //std::thread([filename]() {
        unsigned char* screenshotData = new unsigned char[windowWidth * windowHeight * 3];
        glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);

        ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData);

        if (screenshotImg.save(filename, ImageIO::FORMAT_JPEG) == ImageIO::OK)
            cout << "File " << filename << " saved successfully." << endl;
        else cout << "Failed to save file " << filename << '.' << endl;
        //Sleep(67);
        delete[] screenshotData;
    //}).detach();
}

void idleFunc()
{
    static auto lastScreenshotTime = std::chrono::high_resolution_clock::now();
    // animation
    static int fileIndex = 0;
    char fileName[20];

    auto currentTime = std::chrono::high_resolution_clock::now(); // current time
    float elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastScreenshotTime).count(); // Elapsed time in milliseconds

    if (takeScreenshot && elapsed > 670) {
        if (fileIndex < 1000) {
            sprintf(fileName, "screenshot/%03d.jpg", fileIndex);
            saveScreenshot(fileName);
            //Sleep(2000); // sleep for around 1/15 second
            fileIndex++;
            lastScreenshotTime = currentTime; // Update the last screenshot time
        }
    
    }

    
    // Notify GLUT that it should call displayFunc.
    glutPostRedisplay();
}

void reshapeFunc(int w, int h)
{
    glViewport(0, 0, w, h);

    // When the window has been resized, we need to re-set our projection matrix.
    matrix.SetMatrixMode(OpenGLMatrix::Projection);
    matrix.LoadIdentity();
    // You need to be careful about setting the zNear and zFar. 
    // Anything closer than zNear, or further than zFar, will be culled.
    const float zNear = 0.1f;
    const float zFar = 10000.0f;
    const float humanFieldOfView = 60.0f;
    matrix.Perspective(humanFieldOfView, 1.0f * w / h, zNear, zFar);
  
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
}

void loadSpline(char* argv)
{
    FILE* fileSpline = fopen(argv, "r");
    if (fileSpline == NULL)
    {
        printf("Cannot open file %s.\n", argv);
        exit(1);
    }

    // Read the number of spline control points.
    fscanf(fileSpline, "%d\n", &spline.numControlPoints);
    printf("Detected %d control points.\n", spline.numControlPoints);

    // Allocate memory.
    spline.points = (Point*)malloc(spline.numControlPoints * sizeof(Point));
    // Load the control points.
    for (int i = 0; i < spline.numControlPoints; i++)
    {
        if (fscanf(fileSpline, "%lf %lf %lf",
            &spline.points[i].x,
            &spline.points[i].y,
            &spline.points[i].z) != 3)
        {
            printf("Error: incorrect number of control points in file %s.\n", argv);
            exit(1);
        }
    }
}

// Multiply C = A * B, where A is a m x p matrix, and B is a p x n matrix.
// All matrices A, B, C must be pre-allocated (say, using malloc or similar).
// The memory storage for C must *not* overlap in memory with either A or B. 
// That is, you **cannot** do C = A * C, or C = C * B. However, A and B can overlap, and so C = A * A is fine, as long as the memory buffer for A is not overlaping in memory with that of C.
// Very important: All matrices are stored in **column-major** format.
// Example. Suppose 
//      [ 1 8 2 ]
//  A = [ 3 5 7 ]
//      [ 0 2 4 ]
//  Then, the storage in memory is
//   1, 3, 0, 8, 5, 2, 2, 7, 4. 
void MultiplyMatrices(int m, int p, int n, const double* A, const double* B, double* C)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double entry = 0.0;
            for (int k = 0; k < p; k++)
                entry += A[k * m + i] * B[j * p + k];
            C[m * j + i] = entry;
        }
    }
}

// calculate the cross product of two vectors (use the Point structure to represent a vector)
Point crossProduct(Point a, Point b) {
    return Point{ a.y * b.z - a.z * b.y,
                  a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x };
}

// normalize a vector (use the Point structure to represent a vector)
Point normalization(Point a) {
    double length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    double scale = 1.0 / length;
    return Point{a.x * scale, a.y * scale, a.z * scale};
}

// add two vectors (use the Point structure to represent a vector)
Point addPoint(Point a, Point b) {
    return Point{ a.x + b.x, a.y + b.y, a.z + b.z };
}

// add three vectors (use the Point structure to represent a vector)
Point add3Point(Point a, Point b, Point c) {
    return Point{ a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z };
}

// multiply a value and a vector (use the Point structure to represent a vector)
Point multiplyPoint(double a, Point b) {
    return Point{ a * b.x, a * b.y, a * b.z };
}

// calculate the normal of a triangle given point a b and c
Point calculateTriangleNormal(Point a, Point b, Point c) {
    Point edge1 = { b.x - a.x, b.y - a.y, b.z - a.z };
    Point edge2 = { c.x - a.x, c.y - a.y, c.z - a.z };
    Point normal = crossProduct(edge1, edge2);
    normal = normalization(normal);
    return normal;
}

// move the camera
void moveCamera() {

    p += cameraSpeed;
    // to avoid an error
    int pp = static_cast<int>(p) % splinePoints.size();
    Point position = splinePoints[static_cast<int>(pp)];
    Point direction = splineNormalizedTangents[static_cast<int>(pp)];
    Point up = splineNormals[static_cast<int>(pp)];
    // calculate the new lookat direction by adding the direction to the position
    Point lookAtPosition = { position.x + direction.x, position.y + direction.y, position.z + direction.z };

    // set the lookat matrix
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.LoadIdentity();
    matrix.LookAt(position.x, position.y, position.z,
                  lookAtPosition.x, lookAtPosition.y, lookAtPosition.z,
                  up.x, up.y, up.z);
}

void mouseMotionDragFunc(int x, int y)
{
    // Mouse has moved, and one of the mouse buttons is pressed (dragging).

    // the change in mouse position since the last invocation of this function
    int mousePosDelta[2] = { x - mousePos[0], y - mousePos[1] };

    switch (controlState)
    {
    // translate the terrain
    case TRANSLATE:
        if (leftMouseButton)
        {
        // control x,y translation via the left mouse button
        terrainTranslate[0] += mousePosDelta[0] * 0.01f;
        terrainTranslate[1] -= mousePosDelta[1] * 0.01f;
        }
        if (middleMouseButton)
        {
        // control z translation via the middle mouse button
        terrainTranslate[2] += mousePosDelta[1] * 0.01f;
        }
        break;

    // rotate the terrain
    case ROTATE:
        if (leftMouseButton)
        {
        // control x,y rotation via the left mouse button
        terrainRotate[0] += mousePosDelta[1];
        terrainRotate[1] += mousePosDelta[0];
        }
        if (middleMouseButton)
        {
        // control z rotation via the middle mouse button
        terrainRotate[2] += mousePosDelta[1];
        }
        break;

    // scale the terrain
    case SCALE:
        if (leftMouseButton)
        {
        // control x,y scaling via the left mouse button
        terrainScale[0] *= 1.0f + mousePosDelta[0] * 0.01f;
        terrainScale[1] *= 1.0f - mousePosDelta[1] * 0.01f;
        }
        if (middleMouseButton)
        {
        // control z scaling via the middle mouse button
        terrainScale[2] *= 1.0f - mousePosDelta[1] * 0.01f;
        }
        break;
    }

    // store the new mouse position
    mousePos[0] = x;
    mousePos[1] = y;
}

void mouseMotionFunc(int x, int y)
{
    // Mouse has moved.
    // Store the new mouse position.
    mousePos[0] = x;
    mousePos[1] = y;
}

void mouseButtonFunc(int button, int state, int x, int y)
{
    // A mouse button has has been pressed or depressed.

    // Keep track of the mouse button state, in leftMouseButton, middleMouseButton, rightMouseButton variables.
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
        leftMouseButton = (state == GLUT_DOWN);
    break;

    case GLUT_MIDDLE_BUTTON:
        middleMouseButton = (state == GLUT_DOWN);
    break;

    case GLUT_RIGHT_BUTTON:
        rightMouseButton = (state == GLUT_DOWN);
    break;
    }

    // Keep track of whether CTRL and SHIFT keys are pressed.
    switch (glutGetModifiers())
    {
    case GLUT_ACTIVE_CTRL:
        controlState = TRANSLATE;
    break;

    case GLUT_ACTIVE_SHIFT:
        controlState = SCALE;
    break;

    // If CTRL and SHIFT are not pressed, we are in rotate mode.
    default:
        controlState = ROTATE;
    break;
    }

    // Store the new mouse position.
    mousePos[0] = x;
    mousePos[1] = y;
}

void keyboardFunc(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '1':
        shaderMode = NOTSMOOTH;
    break;

    case 27: // ESC key
        exit(0); // exit the program
    break;

    case 'a':
        takeScreenshot = true;
    break;

    case ' ':
        cout << "You pressed the spacebar." << endl;
    break;

    case 'x':
        // Take a screenshot.
        saveScreenshot("screenshot.jpg");
    break;
    }
}

// to do this: viewLightDirection = (view * float4(lightDirection, 0.0)).xyz;
void transformVector(const float* matrix, const float* vector, float* newVector) {
    float x = vector[0], y = vector[1], z = vector[2], w = 0.0f;
    newVector[0] = matrix[0] * x + matrix[4] * y + matrix[8] * z + matrix[12] * w;
    newVector[1] = matrix[1] * x + matrix[5] * y + matrix[9] * z + matrix[13] * w;
    newVector[2] = matrix[2] * x + matrix[6] * y + matrix[10] * z + matrix[14] * w;
}

void displayFunc()
{
    // This function performs the actual rendering.

    // First, clear the screen.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up the camera position, focus point, and the up vector.
    
    //matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    //matrix.LoadIdentity();
    //matrix.LookAt(6.0, 0.0, 15.0,
    //    10.0, 0.0, 0.0,
    //    0.0, 1.0, 0.0);
    moveCamera();

    // perform translations, rotations and scales.
    matrix.Translate(terrainTranslate[0], terrainTranslate[1], terrainTranslate[2]);
    matrix.Rotate(terrainRotate[0], 1.0, 0.0, 0.0);
    matrix.Rotate(terrainRotate[1], 0.0, 1.0, 0.0);
    matrix.Rotate(terrainRotate[2], 0.0, 0.0, 1.0);
    matrix.Scale(terrainScale[0], terrainScale[1], terrainScale[2]);

    pipelineProgram->Bind();

    // Read the current modelview and projection matrices from our helper class.
    // The matrices are only read here; nothing is actually communicated to OpenGL yet.
    float modelViewMatrix[16];
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.GetMatrix(modelViewMatrix);

    float projectionMatrix[16];
    matrix.SetMatrixMode(OpenGLMatrix::Projection);
    matrix.GetMatrix(projectionMatrix);

    float view[16];
    matrix.GetMatrix(view);

    // get a handle to the program
    GLuint programHandle = pipelineProgram->GetProgramHandle();
    // get a handle to the normalMatrix shader variable
    GLint h_viewLightDirection = glGetUniformLocation(programHandle, "viewLightDirection");

    float lightDirection[3] = { 0.0, -1.0, 0 }; // the “Sun” at noon
    float viewLightDirection[3]; // light direction in the view space
    // the following line is pseudo-code:
    //viewLightDirection = (view * float4(lightDirection, 0.0)).xyz;
    transformVector(view, lightDirection, viewLightDirection);
    // upload viewLightDirection to the GPU
    glUniform3fv(h_viewLightDirection, 1, viewLightDirection);

    // get a handle to the program
    GLuint GetProgramHandle = pipelineProgram->GetProgramHandle();
    // get a handle to the normalMatrix shader variable
    GLint h_normalMatrix = glGetUniformLocation(GetProgramHandle, "normalMatrix");
    float n[16];
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.GetNormalMatrix(n); // get normal matrix
    // upload n to the GPU
    GLboolean isRowMajor = GL_FALSE;
    glUniformMatrix4fv(h_normalMatrix, 1, isRowMajor, n);


    pipelineProgram->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
    pipelineProgram->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
    pipelineProgram->SetUniformVariableMatrix4fv("normalMatrix", GL_FALSE, view);

    vao1->Bind();
    // draw the track tube by tube
    for (int i = 0;i < splinePoints.size() - 1; ++i) {
        glDrawArrays(GL_TRIANGLES, i * numTubeVertices, numTubeVertices);
    }


    // ground
    glActiveTexture(GL_TEXTURE0); 
    
    pipelineProgramGround->Bind();
    glBindTexture(GL_TEXTURE_2D, textureHandle);
    vaoGround->Bind();
    pipelineProgramGround->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
    pipelineProgramGround->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);

    glDrawArrays(GL_TRIANGLES, 0, 6);

    // Swap the double-buffers.
    glutSwapBuffers();
}

vector<Point> calculateSplinePoints(const Spline& spline, int numPointsPerSegment) {
    vector<Point> points;

    const double M[16] = { -0.5, 1.0, -0.5, 0.0,
                           1.5, -2.5, 0.0, 1.0,
                           -1.5, 2.0, 0.5, 0.0,
                           0.5, -0.5, 0.0, 0.0
    };

    for (int i = 0; i <= spline.numControlPoints - 4; i++) {
        double C[12]; // control points
        for (int j = 0; j < 4; j++) {
            C[j + 0 * 4] = spline.points[i + j].x; // column 1, 2, 3, 4 (x)
            C[j + 1 * 4] = spline.points[i + j].y; // column 1, 2, 3, 4 (y)
            C[j + 2 * 4] = spline.points[i + j].z; // column 1, 2, 3, 4 (z)
        }

        for (float u = 0; u <= 1.0; u += 1.0 / numPointsPerSegment) {
            double U[4] = { u * u * u, u * u, u, 1.0 };
            double UmM[4], point[3];
            MultiplyMatrices(1, 4, 4, U, M, UmM);
            MultiplyMatrices(1, 4, 3, UmM, C, point);
            points.push_back({ point[0], point[1], point[2]});
        }
    }
    return points;
}

vector<Point> calculateSplineTangents(const Spline& spline, int numPointsPerSegment) {
    vector<Point> tangents;
    const double M[16] = { -0.5, 1.0, -0.5, 0.0,
                           1.5, -2.5, 0.0, 1.0,
                           -1.5, 2.0, 0.5, 0.0,
                           0.5, -0.5, 0.0, 0.0 };

    for (int i = 0; i <= spline.numControlPoints - 4; i++) {
        double C[12]; // 
        for (int j = 0; j < 4; j++) {
            C[j + 0 * 4] = spline.points[i + j].x; // x
            C[j + 1 * 4] = spline.points[i + j].y; // y
            C[j + 2 * 4] = spline.points[i + j].z; // z
        }

        for (float u = 0; u <= 1.0; u += 1.0 / numPointsPerSegment) {
            double U[4] = { 3 * u * u, 2 * u, 1, 0 }; // Derivative of U vector for tangents
            double UmM[4], tangent[3];
            MultiplyMatrices(1, 4, 4, U, M, UmM); // U'M
            MultiplyMatrices(1, 4, 3, UmM, C, tangent); // U'MC
            tangents.push_back({ tangent[0], tangent[1], tangent[2] });
        }
    }
    return tangents;
}

void calculateNormalsAndBinormals(vector<Point>& tangents, vector<Point>& normals, vector<Point>& binormals) {
    Point arbitraryV = { 0, 0, -1 }; // an arbitrary vector

    normals.resize(tangents.size());
    binormals.resize(tangents.size());

    // initial normal and binormal vectors
    Point T0 = tangents[0];
    Point N0 = normalization(crossProduct(T0, arbitraryV));
    Point B0 = normalization(crossProduct(T0, N0));
    normals[0] = N0;
    binormals[0] = B0;

    // calculate N and B for each point on the spline
    for (int i = 1; i < tangents.size(); i++) {
        Point T = tangents[i];
        Point N = normalization(crossProduct(B0, T));
        Point B = normalization(crossProduct(T, N));

        normals[i] = N;
        binormals[i] = B;

        B0 = B;
    }
}

int initTexture(const char* imageFilename, GLuint textureHandle)
{
    // Read the texture image.
    ImageIO img;
    ImageIO::fileFormatType imgFormat;
    ImageIO::errorType err = img.load(imageFilename, &imgFormat);
    if (err != ImageIO::OK)
    {
        printf("Loading texture from %s failed.\n", imageFilename);
        return -1;
    }
    // Check that the number of bytes is a multiple of 4.
    if (img.getWidth() * img.getBytesPerPixel() % 4)
    {
        printf("Error (%s): The width*numChannels in the loaded image must be a multiple of 4.\n", imageFilename);
        return -1;
    }
    // Allocate space for an array of pixels.
    int width = img.getWidth();
    int height = img.getHeight();
    unsigned char* pixelsRGBA = new unsigned char[4 * width * height]; // we will use 4 bytes per pixel, i.e., RGBA

    // Fill the pixelsRGBA array with the image pixels.
    memset(pixelsRGBA, 0, 4 * width * height); // set all bytes to 0
    for (int h = 0; h < height; h++)
        for (int w = 0; w < width; w++)
        {
            // assign some default byte values (for the case where img.getBytesPerPixel() < 4)
            pixelsRGBA[4 * (h * width + w) + 0] = 0; // red
            pixelsRGBA[4 * (h * width + w) + 1] = 0; // green
            pixelsRGBA[4 * (h * width + w) + 2] = 0; // blue
            pixelsRGBA[4 * (h * width + w) + 3] = 255; // alpha channel; fully opaque

            // set the RGBA channels, based on the loaded image
            int numChannels = img.getBytesPerPixel();
            for (int c = 0; c < numChannels; c++) // only set as many channels as are available in the loaded image; the rest get the default value
                pixelsRGBA[4 * (h * width + w) + c] = img.getPixel(w, h, c);
        }

    // Bind the texture.
    glBindTexture(GL_TEXTURE_2D, textureHandle);

    // Initialize the texture.
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelsRGBA);

    // Generate the mipmaps for this texture.
    glGenerateMipmap(GL_TEXTURE_2D);

    // Set the texture parameters.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // Query support for anisotropic texture filtering.
    GLfloat fLargest;
    glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest);
    printf("Max available anisotropic samples: %f\n", fLargest);
    // Set anisotropic texture filtering.
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 0.5f * fLargest);

    // Query for any errors.
    GLenum errCode = glGetError();
    if (errCode != 0)
    {
        printf("Texture initialization error. Error code: %d.\n", errCode);
        return -1;
    }

    // De-allocate the pixel array -- it is no longer needed.
    delete[] pixelsRGBA;

    return 0;
}

void setPositionandColor(Point v, Point normal, int indexOfVertix) {
    positions[indexOfVertix * 3 + 0] = v.x;
    positions[indexOfVertix * 3 + 1] = v.y;
    positions[indexOfVertix * 3 + 2] = v.z;
    // the varient name is wierd because there are a lot of normals here, this is to avoid getting confused by similar names
    colortonormals[indexOfVertix * 3 + 0] = normal.x; // Red
    colortonormals[indexOfVertix * 3 + 1] = normal.y; // Green
    colortonormals[indexOfVertix * 3 + 2] = normal.z; // Blue
    //normals[indexOfVertix * 4 + 3] = 1.0f; //
    indexOfVertix++;
}

void initScene(int argc, char *argv[])
{
    
    // Set the background color.
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black color.

    // Enable z-buffering (i.e., hidden surface removal using the z-buffer algorithm).
    glEnable(GL_DEPTH_TEST);


    // TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK--TRACK


    // Create a pipeline program. This operation must be performed BEFORE we initialize any VAOs.
    pipelineProgram = new PipelineProgram(); // Load and set up the pipeline program, including its shaders.
    // Load and set up the pipeline program, including its shaders.
    if (pipelineProgram->BuildShadersFromFiles(shaderBasePath, "vertexShader.glsl", "fragmentShader.glsl") != 0)
    {
        cout << "Failed to build the pipeline program." << endl;
        throw 1;
    } 
    cout << "Successfully built the pipeline program." << endl;
    
    // bind the pipeline program
    pipelineProgram->Bind();

    int indexOfVertix = 0;
    // calculate spline points
    splinePoints = calculateSplinePoints(spline, 100); // 100 points per segment

    // calculate tangents, normals and binormals
    splineTangents = calculateSplineTangents(spline, 100);
    calculateNormalsAndBinormals(splineTangents, splineNormals, splineBinormals);

    // normalize tanget
    for (const auto& splineTangent : splineTangents) {
        splineNormalizedTangents.push_back(normalization(splineTangent));
    }

    numTubeVertices = 48;
    numTube = splinePoints.size()-1;
    numVertices = numTube * numTubeVertices*6;
    // vertex positions.
    positions = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex
    // vertex colors.
    colortonormals = (float*)malloc(numVertices * 3 * sizeof(float)); // 4 floats per vertex

    // set positions and colors

    for (size_t i = 0; i < splinePoints.size() -1; ++i) {

        Point sp = splinePoints[i];
        Point sn = splineNormals[i];
        Point sb = splineBinormals[i]; 

        Point sp1 = splinePoints[i+1];
        Point sn1 = splineNormals[i+1];
        Point sb1 = splineBinormals[i+1];

        // see the shape of cross section in the readme
        Point v0 = add3Point(sp, multiplyPoint(-1.0 * unit, sn), multiplyPoint(-2.0 * unit, sb));
        Point v1 = add3Point(sp, multiplyPoint(-1.0 * unit, sn), multiplyPoint(-1.0 * unit, sb));
        Point v2 = add3Point(sp, multiplyPoint(-2.0 * unit, sn), multiplyPoint(-1.0 * unit, sb));
        Point v3 = add3Point(sp, multiplyPoint(-2.0 * unit, sn), multiplyPoint(1.0 * unit, sb));
        Point v4 = add3Point(sp, multiplyPoint(-1.0 * unit, sn), multiplyPoint(1.0 * unit, sb));
        Point v5 = add3Point(sp, multiplyPoint(-1.0 * unit, sn), multiplyPoint(2.0 * unit, sb));
        Point v6 = add3Point(sp, multiplyPoint(-3.0 * unit, sn), multiplyPoint(2.0 * unit, sb));
        Point v7 = add3Point(sp, multiplyPoint(-3.0 * unit, sn), multiplyPoint(-2.0 * unit, sb));

        Point v8 = add3Point(sp1, multiplyPoint(-1.0 * unit, sn1), multiplyPoint(-2.0 * unit, sb1));
        Point v9 = add3Point(sp1, multiplyPoint(-1.0 * unit, sn1), multiplyPoint(-1.0 * unit, sb1));
        Point v10 = add3Point(sp1, multiplyPoint(-2.0 * unit, sn1), multiplyPoint(-1.0 * unit, sb1));
        Point v11 = add3Point(sp1, multiplyPoint(-2.0 * unit, sn1), multiplyPoint(1.0 * unit, sb1));
        Point v12 = add3Point(sp1, multiplyPoint(-1.0 * unit, sn1), multiplyPoint(1.0 * unit, sb1));
        Point v13 = add3Point(sp1, multiplyPoint(-1.0 * unit, sn1), multiplyPoint(2.0 * unit, sb1));
        Point v14 = add3Point(sp1, multiplyPoint(-3.0 * unit, sn1), multiplyPoint(2.0 * unit, sb1));
        Point v15 = add3Point(sp1, multiplyPoint(-3.0 * unit, sn1), multiplyPoint(-2.0 * unit, sb1));


        /* node index in vbo, draw triangle. it would lead to some problem if draw triangle strip
        0 8 1, 8 1 9, 1 9 2, 9 2 10, 2 10 3, 10 3 11, 3 11 4, 11 4 12,
        4 12 5, 12 5 13, 5 13 6, 13 6 14, 6 14 7, 14 7 15, 7 15 0, 15 0 8
        */

        //1
        setPositionandColor(v0, calculateTriangleNormal(v0, v8, v1), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v8, calculateTriangleNormal(v0, v8, v1), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v1, calculateTriangleNormal(v0, v8, v1), indexOfVertix);
        indexOfVertix++;
        //2
        setPositionandColor(v8, calculateTriangleNormal(v0, v8, v1), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v1, calculateTriangleNormal(v0, v8, v1), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v9, calculateTriangleNormal(v0, v8, v1), indexOfVertix);
        indexOfVertix++;
        //3
        setPositionandColor(v1, calculateTriangleNormal(v1, v9, v2), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v9, calculateTriangleNormal(v1, v9, v2), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v2, calculateTriangleNormal(v1, v9, v2), indexOfVertix);
        indexOfVertix++;
        //4
        setPositionandColor(v9, calculateTriangleNormal(v1, v9, v2), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v2, calculateTriangleNormal(v1, v9, v2), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v10, calculateTriangleNormal(v1, v9, v2), indexOfVertix);
        indexOfVertix++;
        //5
        setPositionandColor(v2, calculateTriangleNormal(v2, v10, v3), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v10, calculateTriangleNormal(v2, v10, v3), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v3, calculateTriangleNormal(v2, v10, v3), indexOfVertix);
        indexOfVertix++;
        //6
        setPositionandColor(v10, calculateTriangleNormal(v2, v10, v3), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v3, calculateTriangleNormal(v2, v10, v3), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v11, calculateTriangleNormal(v2, v10, v3), indexOfVertix);
        indexOfVertix++;
        //7
        setPositionandColor(v3, calculateTriangleNormal(v3, v11, v4), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v11, calculateTriangleNormal(v3, v11, v4), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v4, calculateTriangleNormal(v3, v11, v4), indexOfVertix);
        indexOfVertix++;
        //8
        setPositionandColor(v11, calculateTriangleNormal(v3, v11, v4), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v4, calculateTriangleNormal(v3, v11, v4), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v12, calculateTriangleNormal(v3, v11, v4), indexOfVertix);
        indexOfVertix++;
        //9
        setPositionandColor(v4, calculateTriangleNormal(v4, v12, v5), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v12, calculateTriangleNormal(v4, v12, v5), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v5, calculateTriangleNormal(v4, v12, v5), indexOfVertix);
        indexOfVertix++;
        //10
        setPositionandColor(v12, calculateTriangleNormal(v4, v12, v5), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v5, calculateTriangleNormal(v4, v12, v5), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v13, calculateTriangleNormal(v4, v12, v5), indexOfVertix);
        indexOfVertix++;
        //11
        setPositionandColor(v5, calculateTriangleNormal(v5, v13, v6), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v13, calculateTriangleNormal(v5, v13, v6), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v6, calculateTriangleNormal(v5, v13, v6), indexOfVertix);
        indexOfVertix++;
        //12
        setPositionandColor(v13, calculateTriangleNormal(v5, v13, v6), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v6, calculateTriangleNormal(v5, v13, v6), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v14, calculateTriangleNormal(v5, v13, v6), indexOfVertix);
        indexOfVertix++;
        //13
        setPositionandColor(v6, calculateTriangleNormal(v6, v14, v7), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v14, calculateTriangleNormal(v6, v14, v7), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v7, calculateTriangleNormal(v6, v14, v7), indexOfVertix);
        indexOfVertix++;
        //14
        setPositionandColor(v14, calculateTriangleNormal(v6, v14, v7), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v7, calculateTriangleNormal(v6, v14, v7), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v15, calculateTriangleNormal(v6, v14, v7), indexOfVertix);
        indexOfVertix++;
        //15
        setPositionandColor(v7, calculateTriangleNormal(v7, v15, v0), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v15, calculateTriangleNormal(v7, v15, v0), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v0, calculateTriangleNormal(v7, v15, v0), indexOfVertix);
        indexOfVertix++;
        //16
        setPositionandColor(v15, calculateTriangleNormal(v7, v15, v0), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v0, calculateTriangleNormal(v7, v15, v0), indexOfVertix);
        indexOfVertix++;
        setPositionandColor(v8, calculateTriangleNormal(v7, v15, v0), indexOfVertix);
        indexOfVertix++;

    }

    // Initialize VBO for vertices
    vboVertices1 = new VBO(numVertices, 3, positions, GL_STATIC_DRAW);

    // Create the VBOs. 
    //vboVertices1 = new VBO(numVertices, 3, positions, GL_STATIC_DRAW); // 3 values per position
    vboNormals1 = new VBO(numVertices, 3, colortonormals, GL_STATIC_DRAW); // 4 values per color
    // Create the VAOs.
    vao1 = new VAO();
    // Set up the relationship between the "position" shader variable and the VAO.
    vao1->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVertices1, "position");
    // Set up the relationship between the "color" shader variable and the VAO.
    vao1->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboNormals1, "normal");

    free(positions);
    free(colortonormals);


    // GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND--GROUND

    glGenTextures(1, &textureHandle);
    if (initTexture("grass.jpg", textureHandle) != 0)
    {
        std::cerr << "Failed to initialize ground texture." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Create a pipeline program. This operation must be performed BEFORE we initialize any VAOs.
    pipelineProgramGround = new PipelineProgram(); // Load and set up the pipeline program, including its shaders.
    // Load and set up the pipeline program, including its shaders.
    if (pipelineProgramGround->BuildShadersFromFiles(shaderBasePath, "vertexShaderGround.glsl", "fragmentShaderGround.glsl") != 0)
    {
        cout << "Failed to build the ground pipeline program." << endl;
        throw 1;
    }
    cout << "Successfully built the ground pipeline program." << endl;

    int indexOfVertixGround = 0;

    numVerticesGround = 6;
    // vertex positions.
    positionsGround = (float*)malloc(numVerticesGround * 3 * sizeof(float)); // 3 floats per vertex
    // vertex colors.
    colorsGround = (float*)malloc(numVerticesGround * 2 * sizeof(float)); // 2 floats per vertex

    // draw two triangles (one square) to represent the ground
    float positionsGroundData[] = {
        // triangle one
        -100.0f, -50.0f, -100.0f, 
         100.0f, -50.0f, -100.0f, 
        -100.0f, -50.0f,  100.0f, 
        // triangle two
         100.0f, -50.0f, -100.0f, 
         100.0f, -50.0f,  100.0f, 
        -100.0f, -50.0f,  100.0f  
    };

    // trxture coord
    float texCoordsGroundData[] = {
        // triangle one
        0.0f, 0.0f,
        1.0f, 0.0f,
        0.0f, 1.0f,
        // triangle two
        1.0f, 0.0f,
        1.0f, 1.0f,
        0.0f, 1.0f
    };


    // copy the above data into the positionsGround and colorsGround
    memcpy(positionsGround, positionsGroundData, sizeof(positionsGroundData));
    memcpy(colorsGround, texCoordsGroundData, sizeof(texCoordsGroundData));


    // Initialize VBO
    vboVerticesGround = new VBO(numVerticesGround, 3, positionsGround, GL_STATIC_DRAW);
    vboColorsGround = new VBO(numVerticesGround, 2, colorsGround, GL_STATIC_DRAW); // 2 values per color
    // Create the VAOs.
    vaoGround = new VAO();
    // Set up the relationship between the "position" shader variable and the VAO.
    vaoGround->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgramGround, vboVerticesGround, "position");
    // Set up the relationship between the "color" shader variable and the VAO.
    vaoGround->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgramGround, vboColorsGround, "texCoord");

    free(positionsGround);
    free(colorsGround);
  
    // Check for any OpenGL errors.
    std::cout << "GL error status is: " << glGetError() << std::endl;

}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
    cout << "The arguments are incorrect." << endl;
    cout << "usage: ./hw1 <heightmap file>" << endl;
    exit(EXIT_FAILURE);
    }
    if (argc < 2)
    {
        printf("Usage: %s <spline file>\n", argv[0]);
        exit(0);
    }

    cout << "Initializing GLUT..." << endl;
    glutInit(&argc,argv);

    cout << "Initializing OpenGL..." << endl;
    loadSpline("splines/rollerCoaster.sp");

    printf("Loaded spline with %d control point(s).\n", spline.numControlPoints);


    #ifdef __APPLE__
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
    #else
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
    #endif

    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);  
    glutCreateWindow(windowTitle);

    cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
    cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
    cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

    #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(windowWidth - 1, windowHeight - 1);
    #endif

    // Tells GLUT to use a particular display function to redraw.
    glutDisplayFunc(displayFunc);
    // Perform animation inside idleFunc.
    glutIdleFunc(idleFunc);
    // callback for mouse drags
    glutMotionFunc(mouseMotionDragFunc);
    // callback for idle mouse movement
    glutPassiveMotionFunc(mouseMotionFunc);
    // callback for mouse button changes
    glutMouseFunc(mouseButtonFunc);
    // callback for resizing the window
    glutReshapeFunc(reshapeFunc);
    // callback for pressing the keys on the keyboard
    glutKeyboardFunc(keyboardFunc);

    // init glew
    #ifdef __APPLE__
    // nothing is needed on Apple
    #else
    // Windows, Linux
    GLint result = glewInit();
    if (result != GLEW_OK)
    {
        cout << "error: " << glewGetErrorString(result) << endl;
        exit(EXIT_FAILURE);
    }
    #endif

    // Perform the initialization.
    initScene(argc, argv);

    // Sink forever into the GLUT loop.
    glutMainLoop();
}

