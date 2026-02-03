/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Ruyi Bao
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#include <cmath>
#include <iostream>
#include<vector>;

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera, in degrees.
#define fov 60.0

double fovInRadians = fov * 3.1415926 / 180.0;
double aspectRatio = 1.0 * (double)WIDTH / (double)HEIGHT;
// Buffer to store the image when saving it to a JPEG.
unsigned char buffer[HEIGHT][WIDTH][3];

// my resolution
int myWidth = WIDTH / 1;
int myHeight = HEIGHT / 1;

double positive = 0.0001;

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct Ray
{
    double origin[3];
    double direction[3];
};

// for convenient calculation
struct Vector3
{
    double a;
    double b;
    double c;
    Vector3 operator-(const Vector3 m) const {
        return { a - m.a, b - m.b, c - m.c };
    }
    Vector3 operator+(const Vector3 m) const {
        return { a + m.a, b + m.b, c + m.c };
    }
    Vector3 operator*(const double m) const {
        return { a * m, b * m, c * m };
    }
    Vector3 operator/(const double m) const {
        return { a / m, b / m, c / m };
    }
};


struct Intersection {
    bool exist; // whether the intersection exist or not
    Vector3 position;
    Vector3 normal;
    double distance;
    bool isInShadow[MAX_LIGHTS]; // if the intersection exist, whether it is in shadow or not

    double color_diffuse[3];
    double color_specular[3];
    double shininess;

    Vector3 color;

    // constructor
    Intersection()
        : exist(false), position({ 0.0, 0.0, 0.0 }), normal({ 0.0, 0.0, 0.0 }), distance(std::numeric_limits<double>::infinity()), shininess(0.0),
        color({ 0.0, 0.0, 0.0 }) {}
};

// used for checking whether the point is in a 2D triangle
struct Vector2
{
    double a;
    double b;
    Vector2 operator-(const Vector3 m) const {
        return { a - m.a, b - m.b };
    }
};

// variables to store rays and primary intersections
Ray* rays = new Ray[WIDTH * HEIGHT];
Intersection* primaryIntersections = new Intersection[WIDTH * HEIGHT];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

double min_positive = 0.00001;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)// object is triangle
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0) //object is sphere
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

// send out primary rays, store them in Ray* rays
void createPrimaryRay() {
    double tanfov2 = tan(fovInRadians / 2.0);
    for (int w = 0; w < myWidth; w++) {
        for (int h = 0; h < myHeight; h++) {
            // send out rays from the origin with the direction to each pixel center
            double directionX = aspectRatio * tanfov2 * ((w + 0.5) * 2 / double(myWidth) - 1);
            double directionY = tanfov2* ((h + 0.5) * 2 / double(myHeight) - 1);
            Ray& oneRay = rays[w * myHeight + h];
            oneRay.origin[0] = 0.0;
            oneRay.origin[1] = 0.0;
            oneRay.origin[2] = 0.0;
            double length = sqrt(directionX * directionX + directionY * directionY + 1 * 1);
            oneRay.direction[0] = directionX / length;
            oneRay.direction[1] = directionY / length;
            oneRay.direction[2] = -1 / length;
        }
    }
}

double dotProduct(const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double dotProductVector3(const Vector3 m, const Vector3 n) {
    return m.a * n.a + m.b * n.b + m.c * n.c;
}

Vector3 crossProduct(const Vector3& m, const Vector3& n) {
    return {
        m.b * n.c - m.c * n.b,  
        m.c * n.a - m.a * n.c,  
        m.a * n.b - m.b * n.a   
    };
}

void vector3ToDouble(const Vector3& v, double* db) {
    db[0] = v.a;
    db[1] = v.b;
    db[2] = v.c;
}

Vector3 doubleToVector3(const double* a) {
    return { a[0], a[1], a[2] };
}

// calculate the intersection of one ray and one sphere
bool intersectSphere(const Ray& ray, const Sphere& sphere, Intersection& intersection) {
    double temp[3] = { ray.origin[0] - sphere.position[0], ray.origin[1] - sphere.position[1], ray.origin[2] - sphere.position[2] };
    double a = 1.0;     
    double b = 2 * (dotProduct(ray.direction, temp));
    double c = dotProduct(temp, temp) - sphere.radius * sphere.radius;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return false;
    }
    else {
        double discriminantSqrt = sqrt(discriminant);
        double t0 = (-b + discriminantSqrt) / 2.0;
        double t1 = (-b - discriminantSqrt) / 2.0;
        double t = 0.0;
        if (t0 > 0 && t1 > 0) { // two positive t
            if (t0 > t1) {
                t = t1;
            }
            else {
                t = t0;
            }
        }
        else if (t0 > 0) {
            t = t0;
        }
        else if (t1 > 0) {
            t = t1;
        }
        else { // two negative t
            return false; 
        }
        intersection.exist = true;
        intersection.distance = t;
        intersection.position = { ray.origin[0] + t * ray.direction[0], ray.origin[1] + t * ray.direction[1], ray.origin[2] + t * ray.direction[2] };
        Vector3 normal = { intersection.position.a - sphere.position[0], intersection.position.b - sphere.position[1], intersection.position.c - sphere.position[2] };
        double length;
        if (sqrt(dotProductVector3(normal, normal)) > 0) { 
            length = sqrt(dotProductVector3(normal, normal));
        }
        intersection.normal = { normal.a / length, normal.b / length, normal.c / length };
        intersection.color_diffuse[0] = sphere.color_diffuse[0];
        intersection.color_diffuse[1] = sphere.color_diffuse[1];
        intersection.color_diffuse[2] = sphere.color_diffuse[2];
        intersection.color_specular[0] = sphere.color_specular[0];
        intersection.color_specular[1] = sphere.color_specular[1];
        intersection.color_specular[2] = sphere.color_specular[2];
        intersection.shininess = sphere.shininess;
        //std::cout << "Got a primary intersection for a sphere" << std::endl;
        return true;
    }
}

double triangleArea(const Vector2& A, const Vector2& B, const Vector2& C) {
    return ((B.a - A.a) * (C.b - A.b) - ((C.a - A.a) * (B.b - A.b)));
}

bool pointInTriangle(const Vector2& V0, const Vector2& V1, const Vector2& V2, const Vector2& P) {
    double ABC = triangleArea(V0, V1, V2);
    double alpha = triangleArea(P, V1, V2) / ABC;
    double beta = triangleArea(V0, P, V2) / ABC;
    double gamma = 1 - alpha - beta;
    if ( alpha >= 0 && beta >= 0 && gamma >= 0) {
        return true;
    }

    return false;
}

double absolute(double a) {
    if (a < 0) {
        return -a;
    }
    return a;
}

// calculate the intersection of one ray and one triangle
bool intersectTriangle(const Ray& ray, const Triangle& triangle, Intersection& intersection) {
    Vector3 v0 = { triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2] };
    Vector3 v1 = { triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2] };
    Vector3 v2 = { triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2] };

    Vector3 normalOfTriangle = crossProduct(v2 - v1, v0 - v1);

    double length = sqrt(normalOfTriangle.a * normalOfTriangle.a + normalOfTriangle.b * normalOfTriangle.b + normalOfTriangle.c * normalOfTriangle.c);
    normalOfTriangle.a = normalOfTriangle.a / length;
    normalOfTriangle.b = normalOfTriangle.b / length;
    normalOfTriangle.c = normalOfTriangle.c / length;

    double absNormalA = absolute(normalOfTriangle.a);
    double absNormalB = absolute(normalOfTriangle.b);
    double absNormalC = absolute(normalOfTriangle.c);

    // check denominator != 0
    double denominator = dotProductVector3(normalOfTriangle, doubleToVector3(ray.direction));
    if (denominator < min_positive && denominator > -min_positive) { // demoninator == 0
        return false;
    }
    
    double numerator = -(dotProductVector3(normalOfTriangle, doubleToVector3(ray.origin))-dotProductVector3(v0, normalOfTriangle));
    double t = (numerator) / (denominator);

    // intersection is behind the image frame
    if (t < 0) {
        return false;
    }

    //std::cout << "Got a primary intersection for a triangle" << std::endl;
    Vector3 P = doubleToVector3(ray.origin) + doubleToVector3(ray.direction) * t;
    // test if the intersection is inside the triangle
    double alpha;
    double beta;
    double gamma;
    double ABC;
    if (absNormalA >= absNormalB && absNormalA >= absNormalC) {
        // project them to yz-plane
        Vector2 Pyz = { P.b, P.c };
        Vector2 v0yz = { v0.b, v0.c };
        Vector2 v1yz = { v1.b, v1.c };
        Vector2 v2yz = { v2.b, v2.c };
        ABC = triangleArea(v0yz, v1yz, v2yz);
        alpha = triangleArea(Pyz, v1yz, v2yz) / ABC;
        beta = triangleArea(v0yz, Pyz, v2yz) / ABC;
        gamma = 1 - alpha - beta;
    }
    else if (absNormalB >= absNormalC && absNormalB >= absNormalA) {
        // project them to xz-plane
        Vector2 Pxz = { P.a, P.c };
        Vector2 v0xz = { v0.a, v0.c };
        Vector2 v1xz = { v1.a, v1.c };
        Vector2 v2xz = { v2.a, v2.c };
        ABC = triangleArea(v0xz, v1xz, v2xz);
        alpha = triangleArea(Pxz, v1xz, v2xz) / ABC;
        beta = triangleArea(v0xz, Pxz, v2xz) / ABC;
        gamma = 1 - alpha - beta;
    }
    else {

        // project them to xy-plane
        Vector2 Pxy = { P.a, P.b };
        Vector2 v0xy = { v0.a, v0.b };
        Vector2 v1xy = { v1.a, v1.b };
        Vector2 v2xy = { v2.a, v2.b };
        ABC = triangleArea(v0xy, v1xy, v2xy);
        alpha = triangleArea(Pxy, v1xy, v2xy) / ABC;
        beta = triangleArea(v0xy, Pxy, v2xy) / ABC;
        gamma = 1 - alpha - beta;
    }



    if (alpha >= 0 && beta >= 0 && gamma >= 0) {

        intersection.exist = true;
        intersection.distance = t;
        intersection.position = P;
        //std::cout << "---------------------" << std::endl;
        //std::cout << alpha << std::endl;
        //std::cout << beta << std::endl;
        //std::cout << gamma << std::endl;
        // barycentric coordinates for the interpolation of normals, diffuse , specular and shininess
        // normal
        Vector3 normalA = { triangle.v[0].normal[0], triangle.v[0].normal[1], triangle.v[0].normal[2] };
        Vector3 normalB = { triangle.v[1].normal[0], triangle.v[1].normal[1], triangle.v[1].normal[2] };
        Vector3 normalC = { triangle.v[2].normal[0], triangle.v[2].normal[1], triangle.v[2].normal[2] };

        Vector3 interpolatedNormal = (normalA * alpha) + (normalB * beta) + (normalC * gamma);
        double normLength = sqrt(dotProductVector3(interpolatedNormal, interpolatedNormal));
        intersection.normal = { interpolatedNormal.a / normLength, interpolatedNormal.b / normLength, interpolatedNormal.c / normLength };

        // diffuse
        Vector3 diffuseA = { triangle.v[0].color_diffuse[0], triangle.v[0].color_diffuse[1], triangle.v[0].color_diffuse[2] };
        Vector3 diffuseB = { triangle.v[1].color_diffuse[0], triangle.v[1].color_diffuse[1], triangle.v[1].color_diffuse[2] };
        Vector3 diffuseC = { triangle.v[2].color_diffuse[0], triangle.v[2].color_diffuse[1], triangle.v[2].color_diffuse[2] };

        Vector3 interpolatedDiffuse = (diffuseA * alpha) + (diffuseB * beta) + (diffuseC * gamma);
        intersection.color_diffuse[0] = interpolatedDiffuse.a;
        intersection.color_diffuse[1] = interpolatedDiffuse.b;
        intersection.color_diffuse[2] = interpolatedDiffuse.c;

        // specular
        Vector3 specularA = { triangle.v[0].color_specular[0], triangle.v[0].color_specular[1], triangle.v[0].color_specular[2] };
        Vector3 specularB = { triangle.v[1].color_specular[0], triangle.v[1].color_specular[1], triangle.v[1].color_specular[2] };
        Vector3 specularC = { triangle.v[2].color_specular[0], triangle.v[2].color_specular[1], triangle.v[2].color_specular[2] };

        Vector3 interpolatedSpecular = (specularA * alpha) + (specularB * beta) + (specularC * gamma);
        intersection.color_specular[0] = interpolatedSpecular.a;
        intersection.color_specular[1] = interpolatedSpecular.b;
        intersection.color_specular[2] = interpolatedSpecular.c;
        //std::cout << "Got a primary intersection for a triangle" << std::endl;

        intersection.shininess = triangle.v[0].shininess * alpha + triangle.v[1].shininess * beta + triangle.v[2].shininess * gamma;
        return true;
    }
}

// find the nearest intersecion for one ray
Intersection getNearestIntersection(const Ray& ray, const Sphere* spheres, const Triangle* triangles) {
    Intersection nearestIntersection;
    nearestIntersection.exist = false;
    nearestIntersection.distance = std::numeric_limits<double>::infinity();

    // check intersections with all spheres
    for (int i = 0; i < num_spheres; i++) {
        Intersection intersection;
        if (intersectSphere(ray, spheres[i], intersection)) {
            if (intersection.exist && intersection.distance < nearestIntersection.distance) {
                nearestIntersection = intersection;
            }
        }
    }
    //std::cout << "calculating nearest intersection" << std::endl;

    // check intersections with all triangles
    for (int j = 0; j < num_triangles; j++) {
        Intersection intersection;
        
        if (intersectTriangle(ray, triangles[j], intersection)) {
            if (intersection.exist && intersection.distance < nearestIntersection.distance) {
                nearestIntersection = intersection;
            }
        }
    }

    return nearestIntersection;
}

// check if the intersection is in shadow of one light
bool shadowChecking(const Intersection& intersection, const Light light, const Sphere* spheres, const Triangle* triangles) {
        Vector3 direction = doubleToVector3(light.position) - intersection.position;
        // normalize the direction
        double length = sqrt(direction.a * direction.a + direction.b * direction.b + direction.c * direction.c);
        if (length != 0) {
            direction.a = direction.a / length;
            direction.b = direction.b / length;
            direction.c = direction.c / length;
        }
        Ray shadowRay;

        double smallDistance = 0.001;

        shadowRay.origin[0] = intersection.position.a + intersection.normal.a * smallDistance;
        shadowRay.origin[1] = intersection.position.b + intersection.normal.b * smallDistance;
        shadowRay.origin[2] = intersection.position.c + intersection.normal.c * smallDistance;

        shadowRay.direction[0] = direction.a;
        shadowRay.direction[1] = direction.b;
        shadowRay.direction[2] = direction.c;

        // check for intersections with all objects
        Intersection tempIntersection;
        //std::cout << "shadowChecking" << std::endl;

        for (int i = 0; i < num_spheres; i++) {
            if (intersectSphere(shadowRay, spheres[i], tempIntersection) && tempIntersection.distance < length) {
                return true; 
            }
        }
        for (int j = 0; j < num_triangles; j++) {
            if (intersectTriangle(shadowRay, triangles[j], tempIntersection) && tempIntersection.distance < length) {
                return true; 
            }
        }

    return false; 
}

// calculate the nearest intersection for each ray and do the shadow checking, store the result in primaryIntersections
void calculateIntersectionForEachRay(const Sphere* spheres, const Triangle* triangles, const Ray* rays, Intersection*& primaryIntersections) {
    for (int w = 0; w < myWidth; w++) {
        for (int h = 0; h < myHeight; h++) {
            Ray oneRay = rays[w * myHeight + h];
            Intersection nearestIntersection = getNearestIntersection(oneRay, spheres, triangles);
            if (nearestIntersection.exist) {
                primaryIntersections[w * myHeight + h] = nearestIntersection;

                // chack if the intersection is in shadow for each light 
                for (int q = 0; q < num_lights; q++) {
                    primaryIntersections[w * myHeight + h].isInShadow[q] = shadowChecking(nearestIntersection, lights[q], spheres, triangles);
                }
            }
        }
    }
}

Vector3 reflect(const Vector3& m, const Vector3& normal) {
    double dotProduct = dotProductVector3(m, normal);
    return {
        m.a - 2 * dotProduct * normal.a, m.b - 2 * dotProduct * normal.b, m.c - 2 * dotProduct * normal.c
    };
}


void color(const Light* lights, Intersection* primaryIntersections) {
    Vector3 cameraPostion = { 0, 0, 0 };
    for (int h = 0; h < myHeight; h++) {
        for (int w = 0; w < myWidth; w++) {
            Intersection& intersection = primaryIntersections[w * myHeight + h];

            if (!intersection.exist) {
                // background color: white
                intersection.color.a = 255;
                intersection.color.b = 255;
                intersection.color.c = 255;
                continue;
            }

            Vector3 finalColor = {0, 0, 0};

            // calculate color from each light
            for (int i = 0; i < num_lights; i++) {
                const Light& light = lights[i];
                bool inShadow = intersection.isInShadow[i];

                if (!inShadow) {
                    Vector3 lightDirection = {light.position[0] - intersection.position.a,
                                              light.position[1] - intersection.position.b,
                                              light.position[2] - intersection.position.c};
                    double length;
                    if (dotProductVector3(lightDirection, lightDirection) > 0) {
                        length = sqrt(dotProductVector3(lightDirection, lightDirection));
                    }
                    if (length > 0) {
                        lightDirection = { lightDirection.a / length, lightDirection.b / length, lightDirection.c / length };
                    }
                    

                    // diffuse
                    double diff = max(dotProductVector3(lightDirection, intersection.normal), 0.0);
                    Vector3 diffuse = {intersection.color_diffuse[0] * diff,
                                       intersection.color_diffuse[1] * diff,
                                       intersection.color_diffuse[2] * diff};

                    // specular
                    Vector3 reflectDirection = reflect(lightDirection, intersection.normal);
                    Vector3 directionToCamera = intersection.position - cameraPostion;
                    double dlength;
                    if (dotProductVector3(directionToCamera, directionToCamera) > 0) {
                        dlength = sqrt(dotProductVector3(directionToCamera, directionToCamera));
                    }
                    if (dlength > 0) {
                        directionToCamera = { directionToCamera.a / dlength, directionToCamera.b / dlength, directionToCamera.c / dlength };
                    }
                    double spec = pow(max(dotProductVector3(reflectDirection, directionToCamera), 0.0), intersection.shininess);
                    Vector3 specular = {intersection.color_specular[0] * spec,
                                        intersection.color_specular[1] * spec,
                                        intersection.color_specular[2] * spec};

                    // add light color
                    finalColor.a += light.color[0] * (diffuse.a + specular.a);
                    finalColor.b += light.color[1] * (diffuse.b + specular.b);
                    finalColor.c += light.color[2] * (diffuse.c + specular.c);
                }
            }

            // add ambient light
            finalColor.a += ambient_light[0];
            finalColor.b += ambient_light[1];
            finalColor.c += ambient_light[2];

            intersection.color.a = min(max(int(finalColor.a * 255), 0), 255);
            intersection.color.b = min(max(int(finalColor.b * 255), 0), 255);
            intersection.color.c = min(max(int(finalColor.c * 255), 0), 255);
        }
    }
}


void display()
{

}

void draw_scene()
{
    createPrimaryRay();
    calculateIntersectionForEachRay(spheres, triangles, rays, primaryIntersections);
    color(lights, primaryIntersections);
    for (unsigned int x = 0; x < WIDTH; x++)
    {
        glPointSize(2.0);
        // Do not worry about this usage of OpenGL. This is here just so that we can draw the pixels to the screen,
        // after their R,G,B colors were determined by the ray tracer.
        glBegin(GL_POINTS);
        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            // A simple R,G,B output for testing purposes.
            // Modify these R,G,B colors to the values computed by your ray tracer.

            //unsigned char r = x % 256; // modify
            //unsigned char g = y % 256; // modify
            //unsigned char b = (x + y) % 256; // modify
            
            unsigned char r = primaryIntersections[x * myHeight + y].color.a;
            unsigned char g = primaryIntersections[x * myHeight + y].color.b; 
            unsigned char b = primaryIntersections[x * myHeight + y].color.c; 

            plot_pixel(x, y, r, g, b);
        }
        glEnd();
        glFlush();
    }
    printf("Ray tracing completed.\n");
    fflush(stdout);
}


void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

