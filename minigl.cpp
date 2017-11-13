/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

/**
 * Added Data Structures
 */

// vertex
struct vertex {
  // For now vertext contains rgb floats to represent color.
  vec3 color; // How do you assign values to datatypes of vec3?
  vec4 pos;   // Similarly, how do you assign values of datatypes in vec4?
};

struct triangle {
    vertex vertex_one;
    vertex vertex_two;
    vertex vertex_three;
};

struct quad {
    vertex vertex_one;
    vertex vertex_two;
    vertex vertex_three;
    vertex vertex_four;
};
///////////////// End of Added Data Structures //////////////////////

/**
 * Added Global Variables
 */
vec3 current_color = vec3(0,0,0);                    // current color of vertix
MGLpoly_mode draw_mode;                              // Variable set by user to indicate what shape is being drawn.
vector<vertex> list_of_verticies;                    // List of vertices pusdhed in from vertex data structure.
vector<triangle> list_of_triangles;                  // List of triangles built from list of verticies.
vector<mat> matrix_stack;
///////////////// End of Added Global Variables //////////////////////

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * FUNCTION: void mglReadPixels (MGLsize width, MGLsize height, MGLpixel *data)
 * ----------------------------------------------------------
 * AUTHOR: Rick Boshae
 * Date Modified: 11/12/17
 *-----------------------------------------------------------
 * DESCRIPTION:
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 *-----------------------------------------------------------
 * STATUS:
 * Working: Yes
 *
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
  // cout << "In mglReadPixels\n";  // Debugging
  // We're given the framebuffer and our goal is just to fill the buffer. We def do not want to clear vertex information.

  triangle current_triangle; // Used to store a copy of a triangle from a 'list_of_triangles'
  // Begin by iterating through the list of triangle.
  // cout << "list_of_triangles.size():" << list_of_triangles.size() << endl;  // Debugging
  for (unsigned int i = 0; i < list_of_triangles.size(); i++){
  //  cout << "In for loop\n";  // Debugging
    // calculate the pixel coordinates of the triangle
    // determine whether a pixel in the image is inside the triangle transformed to the pixel coordinates
    // you can transform a vertex (x,y) of a triangle to pixel coordinates using
      // i = (x + 1) * width / 2;
      // j = (y + 1) * height / 2;
    current_triangle = list_of_triangles[i];
    // Transform 'current_triangle' from object space to screen space
      // Vertex One Transform
    current_triangle.vertex_one.pos[0] = (current_triangle.vertex_one.pos[0]+1) * width / 2;
    current_triangle.vertex_one.pos[1] = (current_triangle.vertex_one.pos[1]+1) * height / 2;
    // current_triangle.vertex_one.pos[2] = (current_triangle.vertex_one.pos[2]+1) * width / 2;
      // Vertex Two Transform
    current_triangle.vertex_two.pos[0] = (current_triangle.vertex_two.pos[0]+1) * width / 2;
    current_triangle.vertex_two.pos[1] = (current_triangle.vertex_two.pos[1]+1) * height / 2;
    // current_triangle.vertex_two.pos[2] = (current_triangle.vertex_two.pos[2]+1) * width / 2;
      // Vertex Three Transform
    current_triangle.vertex_three.pos[0] = (current_triangle.vertex_three.pos[0]+1) * width / 2;
    current_triangle.vertex_three.pos[1] = (current_triangle.vertex_three.pos[1]+1) * height / 2;
    // current_triangle.vertex_three.pos[2] = (current_triangle.vertex_three.pos[2]+1) * width / 2;

    // Create the bounding box for 'current_triangle'. To create the bounding box find Xmin, Xmax, Ymin, Ymax
    MGLint bounding_box_x_min = (MGLint)floor(min(current_triangle.vertex_one.pos[0], min(current_triangle.vertex_two.pos[0],current_triangle.vertex_three.pos[0]))); // Find Xmin
    MGLint bounding_box_x_max = (MGLint)ceil(max(current_triangle.vertex_one.pos[0], max(current_triangle.vertex_two.pos[0],current_triangle.vertex_three.pos[0]))); // Find Xmax
    MGLint bounding_box_y_min = (MGLint)floor(min(current_triangle.vertex_one.pos[1], min(current_triangle.vertex_two.pos[1],current_triangle.vertex_three.pos[1]))); // Find Ymin
    MGLint bounding_box_y_max = (MGLint)ceil(max(current_triangle.vertex_one.pos[1], max(current_triangle.vertex_two.pos[1],current_triangle.vertex_three.pos[1]))); // Find Ymax

    // Determine the area of 'current_triangle'. Since I am working in two-dimensions, for now I will store the x and y values from each vertex in a new variable.
    vec3 vector_one_to_vector_two    = vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], current_triangle.vertex_two.pos[3]) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], current_triangle.vertex_one.pos[2]);

    vec3 vector_one_to_vector_three  = vec3(current_triangle.vertex_three.pos[0], current_triangle.vertex_three.pos[1], current_triangle.vertex_three.pos[2]) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], current_triangle.vertex_one.pos[2]);

    vec3 vector_two_to_vector_three  = vec3(current_triangle.vertex_three.pos[0], current_triangle.vertex_three.pos[1], current_triangle.vertex_three.pos[2]) - vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], current_triangle.vertex_two.pos[2]);

    // Compute the area 'current_triangle'
    float area_of_triangle = ((cross(vector_one_to_vector_two,vector_two_to_vector_three)).magnitude())/2.0;
    // cout << "vector_one_to_vector_two.magnitude(): " << vector_one_to_vector_two.magnitude() << endl;      // Debugging
    // cout << "vector_one_to_vector_three.magnitude(): " << vector_one_to_vector_three.magnitude() << endl;  // Debugging
    // cout << "area_of_triangle: " << area_of_triangle << endl;                                              // Debugging
    // Iterate through the bounding box to decide whether to draw the pixel and what color it is.

    // Print Where it should be Debugging
    // *(data + (MGLint)current_triangle.vertex_one.pos[0] + (MGLint)current_triangle.vertex_one.pos[1] * width) = Make_Pixel(255,0,0);     // Debugging
    // *(data + (MGLint)current_triangle.vertex_two.pos[0] + (MGLint)current_triangle.vertex_two.pos[1] * width) = Make_Pixel(255,0,0);     // Debugging
    // *(data + (MGLint)current_triangle.vertex_three.pos[0] + (MGLint)current_triangle.vertex_three.pos[1] * width) = Make_Pixel(255,0,0); // Debugging

    for (MGLint y_point = bounding_box_y_min; y_point < bounding_box_y_max; y_point++){
      for (MGLint x_point = bounding_box_x_min; x_point < bounding_box_x_max; x_point++) {
        //
        // the barycentric coordinates can be calculated as

          vec3 vector_one_to_point = vec3(x_point, y_point, 0) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], 0);
          vec3 vector_two_to_point = vec3(x_point, y_point, 0) - vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], 0);
          // alpha = area(Point,vertex_two,vertex_three) / area(vertex_one,vertex_two,vertex_three)
          float alpha = ((cross(vector_two_to_vector_three, vector_two_to_point)).magnitude())/(2 * area_of_triangle);
          // beta = area(vertex_one,Point,vertex_three) / area(vertex_one,vertex_two,vertex_three)
          float beta  = ((cross(vector_one_to_vector_three, vector_one_to_point)).magnitude())/(2 * area_of_triangle);
          // gamma = area(vertex_one,vertex_two, Point) / area(vertex_one,vertex_two,vertex_three);
          float gamma = ((cross(vector_one_to_vector_two, vector_one_to_point)).magnitude())/(2 * area_of_triangle);

           // cout << "alpha , beta , gamma: " << alpha << " " <<  beta << " " << " " << gamma << endl; // Debugging
           //cout << "alpha + beta + gamma: " << alpha + beta + gamma << endl;                          // Debugging
          if(alpha + beta + gamma <= 1.0) {
            // cout << "LESS THAN 1 ! alpha + beta + gamma: " << alpha + beta + gamma << endl;          // Debugging
            *(data + x_point + y_point * width) = Make_Pixel(255,255,255); // QUESTION: How do I read this part of the code?
          } else {
            //*(data + x_point + y_point * width) = Make_Pixel(100,0,0);                                 // Debugging
          }

      }
    }
  }

} // End of mglReadPixels

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
  draw_mode = mode;          // Set draw mode specified by user.
  //list_of_verticies.clear();   // Be nice and clear list of verticies for user.
  //list_of_triangles.clear();   // Be nice and clear list of traingles for user.
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  switch (draw_mode) {
    // Convert list_of_verticies into list_of_triangles.
    case MGL_TRIANGLES :
      for(unsigned int i = 0; (i + 3) <= list_of_verticies.size(); i+=3) {
        triangle current_triangle;   // current_triangle represents the triangle about the be built
        current_triangle.vertex_one   = list_of_verticies.at(i);
        current_triangle.vertex_two   = list_of_verticies.at(i + 1);
        current_triangle.vertex_three = list_of_verticies.at(i + 2);

        // Push fully built triangle onto list_of_triangles vector.
        list_of_triangles.push_back(current_triangle);
      }
    break;

    // Convert list_of_verticies (which represents Quads) into list_of_triangles.

    // NAIVE IMPLEMENTATION of quads. This way of turning quads to triangles is not ideal. There are cases where the triangles pulled from the quad do not represent the original quad. See example below.
    /*
          1-----4
           \  /  \
            2-----3
    */
    // For instance if we grabbed a,b,c and a,c,d that would be wrong. The program should be made so the user does not have to worry about order of input.

    // UPDATE:
    // The user should expect to start at the first input vertext and draw the next one.
    // The standard convention for describing a shape is from left to right.
    case MGL_QUADS :
      for(unsigned int i = 0; (i + 4) <= list_of_verticies.size(); i+=4) {
        triangle current_triangle;   // first_triangle represents the firt triangle built from the quad.

        // Grab first triangle from quad and store it in current triangle
        current_triangle.vertex_one   = list_of_verticies.at(i);
        current_triangle.vertex_two   = list_of_verticies.at(i + 1);
        current_triangle.vertex_three = list_of_verticies.at(i + 2);

        // Add current_triangle to list of triangles
        list_of_triangles.push_back(current_triangle);

        // Grab first triangle from quad and store it in current triangle
        current_triangle.vertex_one   = list_of_verticies.at(i);
        current_triangle.vertex_two   = list_of_verticies.at(i + 2);
        current_triangle.vertex_three = list_of_verticies.at(i + 3);

        // Add  second current_triangle to list of triangles
        list_of_triangles.push_back(current_triangle);
    }
    break;
  }

  list_of_verticies.clear(); // After 'list_of_verticies' has been converted to 'list_of_triangles'. Clear 'list_of_vertices'.
}



/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
  vertex v;

  v.color = current_color;
  v.pos = vec4(x,y,0,0);   // Still using a four-dimensional vector. z and w are 0.

  list_of_verticies.push_back(v);

}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  vertex v;

  v.color = current_color;
  v.pos = vec4(x,y,z,1);

  list_of_verticies.push_back(v);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{

}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
}
