/*
  CSCI 420 Computer Graphics, Computer Science, USC
  Assignment 1: Height Fields with Shaders.
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

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROL_STATE;
CONTROL_STATE controlState = ROTATE;

typedef enum { ONE_POINT, TWO_LINE, THREE_TRIANGLE, FOUR_SMOOTH } RENDERING_MODE;
RENDERING_MODE renderingMode;

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
char windowTitle[512] = "CSCI 420 Homework 1";

// Stores the image loaded from disk.
ImageIO * heightmapImage;
int heightmapImageHeight = 0;
int heightmapImageWidth = 0;
// Number of vertices in the single triangle (starter code).
int numVertices;

// CSCI 420 helper classes.
OpenGLMatrix matrix;
PipelineProgram * pipelineProgram = nullptr;

// 1
VBO * vboVertices1 = nullptr;
VBO * vboColors1 = nullptr;
VAO * vao1 = nullptr;

// 2
VBO* vboVertices2 = nullptr;
VBO* vboColors2 = nullptr;
VAO* vao2 = nullptr;

// 3
VBO* vboVertices3 = nullptr;
VBO* vboColors3 = nullptr;
VAO* vao3 = nullptr;

// 4
VBO* vboVertices4 = nullptr;
VBO* vboColors4 = nullptr;
VAO* vao4 = nullptr;
VBO* vboVerticesLeft = nullptr;
VBO* vboVerticesRight = nullptr;
VBO* vboVerticesUp = nullptr;
VBO* vboVerticesDown = nullptr;


// Write a screenshot to the specified filename.
void saveScreenshot(const char * filename)
{
    unsigned char * screenshotData = new unsigned char[windowWidth * windowHeight * 3];
    glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);

    ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData);

    if (screenshotImg.save(filename, ImageIO::FORMAT_JPEG) == ImageIO::OK)
        cout << "File " << filename << " saved successfully." << endl;
    else cout << "Failed to save file " << filename << '.' << endl;

    delete [] screenshotData;
}

void idleFunc()
{
    // animation
    static int fileIndex = 0;
    char fileName[8];
    if(fileIndex<300){
        sprintf(fileName, "%03d.jpg", fileIndex);
        saveScreenshot(fileName);
        Sleep(67); // sleep for around 1/15 second
        fileIndex++;
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
        renderingMode = ONE_POINT;
        shaderMode = NOTSMOOTH;
        pipelineProgram->SetUniformVariablei("mode", 0);
        break;
    case '2':
        renderingMode = TWO_LINE;
        shaderMode = NOTSMOOTH;
        pipelineProgram->SetUniformVariablei("mode", 0);
        break;
    case '3':
        renderingMode = THREE_TRIANGLE;
        shaderMode = NOTSMOOTH;
        pipelineProgram->SetUniformVariablei("mode", 0);
        break;
    case '4':
        renderingMode = FOUR_SMOOTH;
        shaderMode = SMOOTH;
        pipelineProgram->SetUniformVariablei("mode", 1);
        break;
    case '+':
        if (shaderMode == SMOOTH) {
            scale *= 2.0;
            pipelineProgram->SetUniformVariablef("scale", scale);
        }
        break;
    case '-':
        if (shaderMode == SMOOTH) {
            scale /= 2.0;
            pipelineProgram->SetUniformVariablef("scale", scale);
        }
        break;
    case '9':
        if (shaderMode == SMOOTH) {
            exponent *= 2.0;
            pipelineProgram->SetUniformVariablef("exponent", exponent);
        }
        break;
    case '0':
        if (shaderMode == SMOOTH) {
            exponent /= 2.0;
            pipelineProgram->SetUniformVariablef("exponent", exponent);
        }
        break;
    case 'j': // turn on jetcolor
        if (shaderMode == SMOOTH) {
            pipelineProgram->SetUniformVariablei("jetColor", 1);
        }
        break;
    case 'k': // turn off jetcolor
        if (shaderMode == SMOOTH) {
            pipelineProgram->SetUniformVariablei("jetColor", 0);
        }
        break;
    case 27: // ESC key
        exit(0); // exit the program
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

void displayFunc()
{
    // This function performs the actual rendering.

    // First, clear the screen.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up the camera position, focus point, and the up vector.
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.LoadIdentity();
    matrix.LookAt(0.5, 1.0, 1.8,
                  0.5, 0.0, 0.0,
                  0.0, 1.0, 0.0); // set camera position

    // perform translations, rotations and scales.
    matrix.Translate(terrainTranslate[0], terrainTranslate[1], terrainTranslate[2]);
    matrix.Rotate(terrainRotate[0], 1.0, 0.0, 0.0);
    matrix.Rotate(terrainRotate[1], 0.0, 1.0, 0.0);
    matrix.Rotate(terrainRotate[2], 0.0, 0.0, 1.0);
    matrix.Scale(terrainScale[0], terrainScale[1], terrainScale[2]);

    // Read the current modelview and projection matrices from our helper class.
    // The matrices are only read here; nothing is actually communicated to OpenGL yet.
    float modelViewMatrix[16];
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.GetMatrix(modelViewMatrix);

    float projectionMatrix[16];
    matrix.SetMatrixMode(OpenGLMatrix::Projection);
    matrix.GetMatrix(projectionMatrix);

    // Upload the modelview and projection matrices to the GPU. Note that these are "uniform" variables.
    // Important: these matrices must be uploaded to *all* pipeline programs used.
    // In hw1, there is only one pipeline program, but in hw2 there will be several of them.
    // In such a case, you must separately upload to *each* pipeline program.
    // Important: do not make a typo in the variable name below; otherwise, the program will malfunction.
    pipelineProgram->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
    pipelineProgram->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);

    // Execute the rendering.
  
    switch (renderingMode) {
    case ONE_POINT:
        vao1->Bind();
        glDrawArrays(GL_POINTS, 0, numVertices);
        break;
    case TWO_LINE:
        vao2->Bind();
        // draw in one direction
        for (int i = 0; i < heightmapImageHeight; i++) {
            glDrawArrays(GL_LINE_STRIP, i * heightmapImageWidth, heightmapImageWidth);
        }
        // draw in the other direction
        for (int j = 0; j < heightmapImageWidth; j++) {
            glDrawArrays(GL_LINE_STRIP, heightmapImageHeight * heightmapImageWidth + j * heightmapImageHeight, heightmapImageHeight);
        }
        break;
    case THREE_TRIANGLE:
        vao3->Bind();
        // Draw the triangle strip row by row
        for (int i = 0; i < heightmapImageHeight-1; i++) {
            glDrawArrays(GL_TRIANGLE_STRIP, i * heightmapImageWidth*2, 2*heightmapImageWidth);
        }
        break;
    case FOUR_SMOOTH:
        vao4->Bind();
        pipelineProgram->SetUniformVariablef("scale", scale);
        pipelineProgram->SetUniformVariablef("exponent", exponent);
        pipelineProgram->SetUniformVariablei("mode", 1);
        // Draw the triangle strip row by row
        for (int i = 0; i < heightmapImageHeight - 1; i++) {
            glDrawArrays(GL_TRIANGLE_STRIP, i * heightmapImageWidth * 2, 2 * heightmapImageWidth);
        }
        break;
    }
  
    // Swap the double-buffers.
    glutSwapBuffers();
}

void initScene(int argc, char *argv[])
{
    // Load the image from a jpeg disk file into main memory.
    heightmapImage = new ImageIO();
    if (heightmapImage->loadJPEG(argv[1]) != ImageIO::OK)
    {
        cout << "Error reading image " << argv[1] << "." << endl;
        exit(EXIT_FAILURE);
    }
    else {
        cout << "Successfully loaded image " << argv[1] << endl;
    }

    // Set the background color.
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black color.

    // Enable z-buffering (i.e., hidden surface removal using the z-buffer algorithm).
    glEnable(GL_DEPTH_TEST);
  
    // Create a pipeline program. This operation must be performed BEFORE we initialize any VAOs.
    pipelineProgram = new PipelineProgram(); // Load and set up the pipeline program, including its shaders.
    // Load and set up the pipeline program, including its shaders.
    if (pipelineProgram->BuildShadersFromFiles(shaderBasePath, "vertexShader.glsl", "fragmentShader.glsl") != 0)
    {
        cout << "Failed to build the pipeline program." << endl;
        throw 1;
    } 
    cout << "Successfully built the pipeline program." << endl;
    
    // Bind the pipeline program that we just created. 
    pipelineProgram->Bind();
  
    heightmapImageHeight = heightmapImage->getHeight();
    heightmapImageWidth = heightmapImage->getWidth();

    // Vertex positions.
    float* positions;

    // Vertex colors.
    float* colors;
    int indexOfVertix = 0;


    // ****************************************************
    // *                                                  *
    // *                    Point                         *
    // *                                                  *
    // ****************************************************

    numVertices = heightmapImageHeight * heightmapImageWidth;

    // Vertex positions.
    positions = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex

    // Vertex colors.
    colors = (float*)malloc(numVertices * 4 * sizeof(float)); // 4 floats per vertex
    indexOfVertix = 0;
    for (int y = 0; y < heightmapImageHeight; y++) {
        for (int x = 0; x < heightmapImageWidth; x++) {
            indexOfVertix = (y * heightmapImageWidth + x);
            positions[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
            positions[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            positions[3 * indexOfVertix + 2] = (1.0) * y / (heightmapImageHeight - 1);
            colors[4 * indexOfVertix] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 2] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 3] = 1.0f;
        }
    }


    // Create the VBOs. 
    vboVertices1 = new VBO(numVertices, 3, positions, GL_STATIC_DRAW); // 3 values per position
    vboColors1 = new VBO(numVertices, 4, colors, GL_STATIC_DRAW); // 4 values per color
    // Create the VAOs.
    vao1 = new VAO();
    // Set up the relationship between the "position" shader variable and the VAO.
    vao1->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVertices1, "position");
    // Set up the relationship between the "color" shader variable and the VAO.
    vao1->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboColors1, "color");
    free(positions);
    free(colors);
    // Check for any OpenGL errors.
    std::cout << "GL error status is: " << glGetError() << std::endl;



    // ****************************************************
    // *                                                  *
    // *                    Line                          *
    // *                                                  *
    // ****************************************************


    numVertices = heightmapImageHeight * heightmapImageWidth*2;
    cout << numVertices << endl;
   
    // Vertex positions.
    positions = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex

    // Vertex colors.
    colors = (float*)malloc(numVertices * 4 * sizeof(float)); // 4 floats per vertex

    // Set all the vertex in one direction
    indexOfVertix = 0;
    for (int y = 0; y < heightmapImageHeight; y++) {
        for (int x = 0; x < heightmapImageWidth; x++) {
            indexOfVertix = (y * heightmapImageWidth + x);
            positions[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
            positions[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            positions[3 * indexOfVertix + 2] = (1.0) * y / (heightmapImageHeight - 1);
            colors[4 * indexOfVertix] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 2] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 3] = 1.0f;
        }
    }
     
    // in the other direction 
    indexOfVertix = heightmapImageHeight * heightmapImageWidth;
    for (int x = 0; x < heightmapImageWidth; x++) {
        for (int y = 0; y < heightmapImageHeight; y++) {
            positions[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
            positions[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            positions[3 * indexOfVertix + 2] = (1.0) * y / (heightmapImageHeight - 1);
            colors[4 * indexOfVertix] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 2] = (1.0) * heightmapImage->getPixel(x, y, 0) / 255.0f;
            colors[4 * indexOfVertix + 3] = 1.0f;
            indexOfVertix += 1;
        }
    }
     
    // Create the VBOs. 
    vboVertices2 = new VBO(numVertices, 3, positions, GL_STATIC_DRAW); // 3 values per position
    vboColors2 = new VBO(numVertices, 4, colors, GL_STATIC_DRAW); // 4 values per color
    // Create the VAOs.
    vao2 = new VAO();
    // Set up the relationship between the "position" shader variable and the VAO.
    vao2->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVertices2, "position");
    // Set up the relationship between the "color" shader variable and the VAO.
    vao2->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboColors2, "color");
    free(positions);
    free(colors);
    // Check for any OpenGL errors.
    std::cout << "GL error status is: " << glGetError() << std::endl;


    // ****************************************************
    // *                                                  *
    // *                    Triangle                      *
    // *                                                  *
    // ****************************************************

    numVertices = heightmapImageHeight * heightmapImageWidth*2; // Actually it's height*width*2-height*2.
    // Vertex positions.
    positions = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex

    // Vertex colors.
    colors = (float*)malloc(numVertices * 4 * sizeof(float)); // 4 floats per vertex

    indexOfVertix = 0;
    for (int y = 0; y < heightmapImageHeight-1; y++) {
        for (int x = 0; x < heightmapImageWidth; x++) {
            for (int m = 0;m <= 1;m++) {
                int newy = y + m;
                positions[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
                positions[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                positions[3 * indexOfVertix + 2] = (1.0) * newy / (heightmapImageHeight - 1);
                colors[4 * indexOfVertix] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                colors[4 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                colors[4 * indexOfVertix + 2] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                colors[4 * indexOfVertix + 3] = 1.0f;
                indexOfVertix += 1;
            }
        }
    }

   
    // Create the VBOs. 
    vboVertices3 = new VBO(numVertices, 3, positions, GL_STATIC_DRAW); // 3 values per position
    vboColors3 = new VBO(numVertices, 4, colors, GL_STATIC_DRAW); // 4 values per color
    // Create the VAOs.
    vao3 = new VAO();
    // Set up the relationship between the "position" shader variable and the VAO.
    vao3->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVertices3, "position");
    // Set up the relationship between the "color" shader variable and the VAO.
    vao3->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboColors3, "color");
    free(positions);
    free(colors);
    // Check for any OpenGL errors.
    std::cout << "GL error status is: " << glGetError() << std::endl;


    // ****************************************************
    // *                                                  *
    // *                    Smooth                        *
    // *                                                  *
    // ****************************************************
    numVertices = heightmapImageHeight * heightmapImageWidth * 2;
    cout << numVertices << endl;
    // Vertex positions.
    positions = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex, i.e., x,y,z
    float* positionsLeft = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex, i.e., x,y,z
    float* positionsRight = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex, i.e., x,y,z
    float* positionsUp = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex, i.e., x,y,z
    float* positionsDown = (float*)malloc(numVertices * 3 * sizeof(float)); // 3 floats per vertex, i.e., x,y,z
    // Vertex colors.
    colors = (float*)malloc(numVertices * 4 * sizeof(float)); // 4 floats per vertex


    indexOfVertix = 0;
    for (int y = 0; y < heightmapImageHeight - 1; y++) {
        for (int x = 0; x < heightmapImageWidth; x++) {
            for (int m = 0;m <= 1;m++) {
                int newy = y + m;
                //centre
                positions[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
                positions[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                positions[3 * indexOfVertix + 2] = (1.0) * newy / (heightmapImageHeight - 1);
                colors[4 * indexOfVertix] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                colors[4 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                colors[4 * indexOfVertix + 2] = (1.0) * heightmapImage->getPixel(x, newy, 0) / 255.0f;
                colors[4 * indexOfVertix + 3] = 1.0f;

                // left  x-1,y
                // right x+1,y
                // up    x,y-1
                // down  x,y+1

                // left  x-1,y
                if (x > 0) {
                    positionsLeft[3 * indexOfVertix] = (1.0) * (x - 1) / (heightmapImageWidth - 1);
                    positionsLeft[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x - 1, newy, 0) / 255.0f;
                    positionsLeft[3 * indexOfVertix + 2] = (1.0) * newy / (heightmapImageHeight - 1);
                }
                else {
                    positionsLeft[3 * indexOfVertix] = positions[3 * indexOfVertix];
                    positionsLeft[3 * indexOfVertix + 1] = positions[3 * indexOfVertix + 1];
                    positionsLeft[3 * indexOfVertix + 2] = positions[3 * indexOfVertix + 2];
                }

                // right x+1,y
                if (x < heightmapImageWidth - 1) {
                    positionsRight[3 * indexOfVertix] = (1.0) * (x + 1) / (heightmapImageWidth - 1);
                    positionsRight[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x + 1, newy, 0) / 255.0f;
                    positionsRight[3 * indexOfVertix + 2] = (1.0) * newy / (heightmapImageHeight - 1);
                }
                else {
                    positionsRight[3 * indexOfVertix] = positions[3 * indexOfVertix];
                    positionsRight[3 * indexOfVertix + 1] = positions[3 * indexOfVertix + 1];
                    positionsRight[3 * indexOfVertix + 2] = positions[3 * indexOfVertix + 2];
                }

                // up    x,y-1
                if (newy > 0) {
                    positionsUp[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
                    positionsUp[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, newy - 1, 0) / 255.0f;
                    positionsUp[3 * indexOfVertix + 2] = (1.0) * (newy - 1) / (heightmapImageHeight - 1);
                }
                else {
                    positionsUp[3 * indexOfVertix] = positions[3 * indexOfVertix];
                    positionsUp[3 * indexOfVertix + 1] = positions[3 * indexOfVertix + 1];
                    positionsUp[3 * indexOfVertix + 2] = positions[3 * indexOfVertix + 2];
                }

                // down  x,y+1
                if (newy < heightmapImageHeight - 1) {
                    positionsDown[3 * indexOfVertix] = (1.0) * x / (heightmapImageWidth - 1);
                    positionsDown[3 * indexOfVertix + 1] = (1.0) * heightmapImage->getPixel(x, newy + 1, 0) / 255.0f;
                    positionsDown[3 * indexOfVertix + 2] = (1.0) * (newy + 1) / (heightmapImageHeight - 1);
                }
                else {
                    positionsDown[3 * indexOfVertix] = positions[3 * indexOfVertix];
                    positionsDown[3 * indexOfVertix + 1] = positions[3 * indexOfVertix + 1];
                    positionsDown[3 * indexOfVertix + 2] = positions[3 * indexOfVertix + 2];
                }

                indexOfVertix += 1;
            }
        }
    }

    // Create the VBOs. 
    vboVertices4 = new VBO(numVertices, 3, positions, GL_STATIC_DRAW); // 3 values per position
    vboVerticesLeft = new VBO(numVertices, 3, positionsLeft, GL_STATIC_DRAW); // 3 values per position
    vboVerticesRight = new VBO(numVertices, 3, positionsRight, GL_STATIC_DRAW); // 3 values per position
    vboVerticesUp = new VBO(numVertices, 3, positionsUp, GL_STATIC_DRAW); // 3 values per position
    vboVerticesDown = new VBO(numVertices, 3, positionsDown, GL_STATIC_DRAW); // 3 values per position
    vboColors4 = new VBO(numVertices, 4, colors, GL_STATIC_DRAW); // 4 values per color
    // Create the VAOs.
    vao4 = new VAO();
    // Set up the relationship between the "position" shader variable and the VAO.
    vao4->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVertices4, "position");
    vao4->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVerticesLeft, "positionLeft");
    vao4->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVerticesRight, "positionRight");
    vao4->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVerticesUp, "positionUp");
    vao4->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVerticesDown, "positionDown");
    // Set up the relationship between the "color" shader variable and the VAO.
    vao4->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboColors4, "color");
    free(positions);
    free(positionsLeft);
    free(positionsRight);
    free(positionsUp);
    free(positionsDown);
    free(colors);
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

    cout << "Initializing GLUT..." << endl;
    glutInit(&argc,argv);

    cout << "Initializing OpenGL..." << endl;

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

