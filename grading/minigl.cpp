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
  vec3 color; // Example of how to assign values to vec3: vec3 color = vec3(0,0,0);
  vec4 pos;   // You can set the value of pos using a simlar example from above.
};

// triangle
struct triangle {
    vertex vertex_one;
    vertex vertex_two;
    vertex vertex_three;
};

// quad
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

vec3 current_color = vec3(0,0,0);                    //'current_color' used to set color of vertix in 'mglVertex2' and 'mglVertex3'
MGLpoly_mode draw_mode;                              // 'draw_mode' is set by user to indicate what shape is being drawn.

vector<vertex> list_of_vertices;                    // 'list_of_vertices' stores vertices in from vertex data structure.
vector<triangle> list_of_triangles;                  // 'list_of_triangles' set from 'list_of_vertices'.


MGLmatrix_mode matrix_mode;                          // 'matrix_mode' is set by user. The matrix mode is either model_view or projection.

mat4 projection_matrix;                              //

mat4 modelview_matrix;                               //

// 'get_current_matrix' is a helper function that returns the corresponding matrix to matrix mode.
mat4& get_current_matrix() {

  switch (matrix_mode) {
    case MGL_PROJECTION:
      return projection_matrix;
    break;

    case MGL_MODELVIEW:
      return modelview_matrix;
    break;

    default :
     cout << "matrix_mode has illegal argument. Exiting program." << endl;
     exit(1);
    break;
  }
}

//mat4 current_matrix;                                 // Current matrix
vector<mat4> model_view_matrix_stack;                  // 'model_view_matrix_stack' stores a list of matricies when matrix mode is set to 'model_view'
vector<mat4> projection_matrix_stack;                  // 'projection_matrix_stack' stores a list of matricies when matrix mode is set to 'projection'

// Function that returns the matrix based on matrix mode.
vector<mat4>& get_current_matrix_stack() {
  if (matrix_mode == MGL_PROJECTION) {
    return projection_matrix_stack;
  } else {
    return model_view_matrix_stack;
  }
}

vector<vector<double>> z_buffer;

///////////////// End of Added Global Variables //////////////////////

/**
 * Helper function
 */
MGLfloat area(vertex a, vertex b, vertex c) {
  return a.pos[0]*(b.pos[1]-c.pos[1]) + a.pos[1]*(c.pos[0]-b.pos[0]) + (b.pos[0]*c.pos[1]-b.pos[1]*c.pos[0]);
}


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
  //cout << "In mglReadPixels\n";  // Debugging

  // Resize z-buffer to match width times height
  z_buffer.resize((MGLint)(width));

  // initialize min-z values in z_buffer
  for (unsigned int buffer_x_index = 0; buffer_x_index < z_buffer.size(); ++buffer_x_index) {
    z_buffer[buffer_x_index].resize(height);

    for (unsigned int buffer_y_index = 0; buffer_y_index < height; ++buffer_y_index )
    z_buffer[buffer_x_index][buffer_y_index] = 2;
  }

  triangle current_triangle; // Used to store a copy of a triangle from a 'list_of_triangles'
//  cout << "list_of_triangles.size() = " << list_of_triangles.size() << endl;
  // Iterate through the 'list_of_triangle'.
  for (unsigned int triangle_list_index = 0; triangle_list_index < list_of_triangles.size(); ++triangle_list_index) {


    current_triangle = list_of_triangles.at(triangle_list_index); // grab triangle from list of triangles.

    // before processing a triangle, coordinates of each vertex is divided by its w to get the normalized device coordinates.
    current_triangle.vertex_one.pos[0] /= current_triangle.vertex_one.pos[3];
    current_triangle.vertex_one.pos[1] /= current_triangle.vertex_one.pos[3];
    current_triangle.vertex_one.pos[2] /= current_triangle.vertex_one.pos[3];

    current_triangle.vertex_two.pos[0] /= current_triangle.vertex_two.pos[3];
    current_triangle.vertex_two.pos[1] /= current_triangle.vertex_two.pos[3];
    current_triangle.vertex_two.pos[2] /= current_triangle.vertex_two.pos[3];

    current_triangle.vertex_three.pos[0] /= current_triangle.vertex_three.pos[3];
    current_triangle.vertex_three.pos[1] /= current_triangle.vertex_three.pos[3];
    current_triangle.vertex_three.pos[2] /= current_triangle.vertex_three.pos[3];

    // cout <<list_of_triangles.at(triangle_list_index).vertex_one.pos[2] << " "<< list_of_triangles.at(triangle_list_index).vertex_two.pos[2]<< " " << list_of_triangles.at(triangle_list_index).vertex_three.pos[2] << endl;
    // cout <<current_triangle.vertex_one.pos[2] << " "<< current_triangle.vertex_two.pos[2]<< " " << current_triangle.vertex_three.pos[2] << endl << endl;
    // calculate the pixel coordinates of the triangle
    // determine whether a pixel in the image is inside the triangle transformed to the pixel coordinates
    // you can transform a vertex (x,y) of a triangle to pixel coordinates using
      // i = (x + 1) * width / 2;
      // j = (y + 1) * height / 2;

    // Vertex One Transform
    current_triangle.vertex_one.pos[0] = (MGLfloat)(((current_triangle.vertex_one.pos[0]+1) * width) / 2); // Transform 'vertex_one' x-coordinate
    current_triangle.vertex_one.pos[1] = (MGLfloat)(current_triangle.vertex_one.pos[1]+1) * height / 2;    // Transform 'vertex_one' y-coordinate

    // Vertex Two Transform
    current_triangle.vertex_two.pos[0] = (MGLfloat)(current_triangle.vertex_two.pos[0]+1) * width / 2;      // Transform 'vertex_two' x-coordinate
    current_triangle.vertex_two.pos[1] = (MGLfloat)(current_triangle.vertex_two.pos[1]+1) * height / 2;     // Transform 'vertex_two' y-coordinate

    // Vertex Three Transform
    current_triangle.vertex_three.pos[0] = (MGLfloat)(current_triangle.vertex_three.pos[0]+1) * width / 2;  // Transform 'vertex_three' x-coordinate
    current_triangle.vertex_three.pos[1] = (MGLfloat)(current_triangle.vertex_three.pos[1]+1) * height / 2; // Transform 'vertex_three' y-coordinate

    // Create the bounding box for 'current_triangle'. To create the bounding box find Xmin, Xmax, Ymin, Ymax
    MGLint bounding_box_x_min = (MGLint)floor(min(current_triangle.vertex_one.pos[0], min(current_triangle.vertex_two.pos[0],current_triangle.vertex_three.pos[0]))); // Find Xmin
    MGLint bounding_box_x_max = (MGLint)ceil(max(current_triangle.vertex_one.pos[0], max(current_triangle.vertex_two.pos[0],current_triangle.vertex_three.pos[0]))); // Find Xmax
    MGLint bounding_box_y_min = (MGLint)floor(min(current_triangle.vertex_one.pos[1], min(current_triangle.vertex_two.pos[1],current_triangle.vertex_three.pos[1]))); // Find Ymin
    MGLint bounding_box_y_max = (MGLint)ceil(max(current_triangle.vertex_one.pos[1], max(current_triangle.vertex_two.pos[1],current_triangle.vertex_three.pos[1]))); // Find Ymax

    if (bounding_box_y_max >= (MGLint) height) {
      bounding_box_y_max = (MGLint) height-1;
    }

    if (bounding_box_y_min < 0) {
      bounding_box_y_min = 0;
    }

    if (bounding_box_x_max >= (MGLint) width) {
      bounding_box_x_max = (MGLint) width-1;
    }

    if (bounding_box_x_min < 0) {
      bounding_box_x_min = 0;
    }
    // Determine the area of 'current_triangle'. Since I am working in two-dimensions, for now I will store the x and y values from each vertex in a new variable.
    // vec3 vertex_one_to_vertex_two    = vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], current_triangle.vertex_two.pos[2]) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], current_triangle.vertex_one.pos[2]);
    //
    // vec3 vertex_one_to_vertex_three  = vec3(current_triangle.vertex_three.pos[0], current_triangle.vertex_three.pos[1], current_triangle.vertex_three.pos[2]) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], current_triangle.vertex_one.pos[2]);
    //
    // vec3 vertex_two_to_vertex_three  = vec3(current_triangle.vertex_three.pos[0], current_triangle.vertex_three.pos[1], current_triangle.vertex_three.pos[2]) - vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], current_triangle.vertex_two.pos[2]);

    // vec3 vertex_one_to_vertex_two    = vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], 0) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], 0);
    //
    // vec3 vertex_one_to_vertex_three  = vec3(current_triangle.vertex_three.pos[0], current_triangle.vertex_three.pos[1], 0) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], 0);
    //
    // vec3 vertex_two_to_vertex_three  = vec3(current_triangle.vertex_three.pos[0], current_triangle.vertex_three.pos[1], 0) - vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], 0);

    // Compute the area 'current_triangle'
    MGLfloat area_of_triangle = area(current_triangle.vertex_one, current_triangle.vertex_two, current_triangle.vertex_three); // I can divide by two here but there is really no need at this point.

    //cout << "area_of_triangle [" << triangle_list_index << "]: " << area_of_triangle/2 << endl << endl;

    for (MGLint y_point   = bounding_box_y_min; y_point <= bounding_box_y_max; ++y_point){
      for (MGLint x_point = bounding_box_x_min; x_point <= bounding_box_x_max; ++x_point) {

        // the barycentric coordinates can be calculated as
          // vec3 vertex_one_to_point = vec3(x_point+.5, y_point+.5, 0) - vec3(current_triangle.vertex_one.pos[0], current_triangle.vertex_one.pos[1], 0);
          // vec3 vertex_two_to_point = vec3(x_point+.5, y_point+.5, 0) - vec3(current_triangle.vertex_two.pos[0], current_triangle.vertex_two.pos[1], 0);

          vertex p;
          p.pos[0] = x_point+.5f;
          p.pos[1] = y_point+.5f;
          //p.pos[2] = 0;           // Not needed
          //p.pos[3] = 1;           // Not needed

          //       alpha = area(Point,vertex_two,vertex_three) / area(vertex_one,vertex_two,vertex_three)
          MGLfloat alpha = area(p, current_triangle.vertex_two, current_triangle.vertex_three)/area_of_triangle;    // I can divide by two here but there is really no need at this point.

          //       beta = area(vertex_one,Point,vertex_three) / area(vertex_one,vertex_two,vertex_three)
          MGLfloat beta  = area(current_triangle.vertex_one, p, current_triangle.vertex_three)/area_of_triangle;    // I can divide by two here but there is really no need at this point.

          //       gamma = area(vertex_one,vertex_two, Point) / area(vertex_one,vertex_two,vertex_three);
          MGLfloat gamma = area(current_triangle.vertex_one, current_triangle.vertex_two, p)/area_of_triangle;      // I can divide by two here but there is really no need at this point.

          //  cout << "alpha , beta , gamma: " << alpha << " " <<  beta << " " << " " << gamma << endl; // Debugging
          //  cout << "alpha + beta + gamma: " << alpha + beta + gamma << endl;                          // Debugging
          if(alpha > 0 && beta > 0 && gamma > 0) {
            // alpha /=area_of_triangle;
            // beta  /=area_of_triangle;
            // gamma /=area_of_triangle;



            MGLfloat alpha_real = (alpha/current_triangle.vertex_one.pos[3])/ ((alpha/current_triangle.vertex_one.pos[3]) + (beta/current_triangle.vertex_two.pos[3]) + (gamma/current_triangle.vertex_three.pos[3]));

            MGLfloat beta_real = (beta/current_triangle.vertex_two.pos[3])/ ((alpha/current_triangle.vertex_one.pos[3]) + (beta/current_triangle.vertex_two.pos[3]) + (gamma/current_triangle.vertex_three.pos[3]));

            MGLfloat gamma_real = (gamma/current_triangle.vertex_three.pos[3])/ ((alpha/current_triangle.vertex_one.pos[3]) + (beta/current_triangle.vertex_two.pos[3]) + (gamma/current_triangle.vertex_three.pos[3]));

            // cout << "alpha:  " << alpha << ", beta: " <<  beta << ", gamma: " << gamma << endl; // Debugging
            // cout << "alpha + beta + gamma: " << alpha + beta + gamma << endl;                          // Debugging

            // FINDING Z DEPTH
            // For finding the depth of a pixel, (i,j), we use the Normalized Device Coordinates (NDC) of
            // our vertices, a,b,c. (in my case 1, 2, 3) I.e, (z′a,z′b,z′c)=(za/wa,zb/wb,zc/wc) and use barycentric coordinates
            // for linear interpolation:
            //
            // Z_{i,j}=α_{i,j}z′_one + β_{i,j}*z′_two + γ_{i,j}z′_three

            // if ((abs(current_triangle.vertex_one.pos[2]) > 1) || (abs(current_triangle.vertex_two.pos[2]) > 1)|| (abs(current_triangle.vertex_three.pos[2]) > 1) ) {
            //   continue;
            // }

            // perform z-interpolation
            //MGLfloat z_depth = alpha_real*(current_triangle.vertex_one.pos[2] * current_triangle.vertex_one.pos[3]) + beta_real*(current_triangle.vertex_two.pos[2]*current_triangle.vertex_two.pos[3]) + gamma_real*(current_triangle.vertex_three.pos[2]*current_triangle.vertex_three.pos[3]);
            double perspective_z_depth = (double)alpha_real*current_triangle.vertex_one.pos[2] + beta_real*current_triangle.vertex_two.pos[2] + gamma_real*current_triangle.vertex_three.pos[2];
            double z_depth = alpha*current_triangle.vertex_one.pos[2] + beta*current_triangle.vertex_two.pos[2] + gamma*current_triangle.vertex_three.pos[2];

            // if((z_depth) > 1) {
            //   //cout << "z-clip too big\n";
            //   //cout << z_depth << endl;
            //   continue;
            // }
            //
            // if(z_depth < -1) {
            //
            //   //cout << "z-clip too big\n";
            //   continue;
            //
            // }

            // cout << "z_depth" << z_depth;

            // check if the z -depth- of the pixel (which we get by
            // z-interpolation) is less than the min_z of that pixel
            // (that we get from our storage). If so, it means that triangle t,
            // at that pixel, is closest and we can color as we did before (set
            // data[i+j*width]) and update the min_z of the pixel.

            //if( (x_point +.5 < 0) || (y_point +.5 < 0) || ( width - (x_point +.5)) <= 0 || ( width - (x_point +.5)) >= width || (height - (y_point+.5)) <= 1 || (height - (y_point+.5)) > height /*|| (z_depth) > 1 || (z_depth < -1)*/)  {
            // if( (x_point < 0) || (y_point < 0) || (y_point >= ((MGLint)height)) || (x_point >= ((MGLint)width)))  {
            //
              // if( (z_depth < 0)) {
              //   cout << "clipped\n";
              //   cout << z_depth << endl;
              // }
            //
            //   continue;
            // }
            //cout << "(height - y_point)" << (height - y_point) << endl;

            if( ((x_point +.5) >= 0) && ((y_point+.5) >= 0) && ((y_point+.5) < ((MGLint)height)) && ((x_point+.5) < ((MGLint)width)) && (z_depth) <= 1.0f && (z_depth >= -1.0f) && (perspective_z_depth <= 1.0f) && (perspective_z_depth >= -1.0f))  {
            // Debuggin -- Seg fault from accessing z-buffer
            // cout << "x_point = " << x_point << ", y_point = " << y_point << ", width = " << width << ", height = " << height << endl;
            // cout  << "z_buffer.size() = " << z_buffer.size() << " Trying to access: " << x_point + y_point * width << endl;
              if(z_buffer[x_point][y_point] > z_depth) {
                // if (z_buffer.at(x_point + y_point * width) != 2) {
                //   cout << "z_depth replaced" << endl;
                //   cout << "old z-value: " << z_buffer.at(x_point + y_point * width) << ", new z-value: " <<  z_depth << endl;
                // }

                // if( z_depth == z_buffer.at(x_point + y_point * width)) {
                //   cout << "z_depth collision" << endl;
                //   cout << "old z-value: " << z_buffer.at(x_point + y_point * width) << ", new z-value: " <<  z_depth << endl;
                // }


                z_buffer[x_point][y_point] = z_depth;


                // Color interpolation
                //vec3 color_interpolation((alpha*current_triangle.vertex_one.color + beta*current_triangle.vertex_two.color + gamma*current_triangle.vertex_three.color)*255);

                // perspective color interpolation
                vec3 color_interpolation((alpha_real*current_triangle.vertex_one.color + beta_real*current_triangle.vertex_two.color + gamma_real*current_triangle.vertex_three.color)*255);

                // Debugging -- Vertex color
                // cout << "Color of Vertex One"  << endl << "R: " << current_triangle.vertex_one.color[0] * 255 << " G: " << current_triangle.vertex_one.color[1] * 255<< " B: " << current_triangle.vertex_one.color[2] * 255 << endl;
                //*(data + x_point + y_point * width) = Make_Pixel(current_triangle.vertex_one.color[0] * 255,current_triangle.vertex_one.color[1] * 255,current_triangle.vertex_one.color[2] * 255 );

                // Debugging -- Vertex color interpolation
                 //cout << "Color interpolation Values"  << endl << "R: " << color_interpolation[0] << " G: " << color_interpolation[1] << " B: " << color_interpolation[2] << endl;
                *(data + x_point + y_point * width) = Make_Pixel(color_interpolation[0], color_interpolation[1],color_interpolation[2]);

              }
            }


      }
    }
  }
}
  //z_buffer.clear(); // I dont think i need this but just in case.
//cout << "Out of mglReadPixels" << endl <<endl; // Debugging
} // End of mglReadPixels

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
  draw_mode = mode;          // Set draw mode specified by user.

}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  //cout << "list_of_vertices.size() " << list_of_vertices.size() << endl;
  switch (draw_mode) {
    // Convert list_of_vertices into list_of_triangles.
    case MGL_TRIANGLES :
      for(unsigned int i = 0; (i + 2) < list_of_vertices.size(); i+=3) {
        triangle current_triangle;   // current_triangle represents the triangle about the be built
        current_triangle.vertex_one   = list_of_vertices.at(i);
        current_triangle.vertex_two   = list_of_vertices.at(i + 1);
        current_triangle.vertex_three = list_of_vertices.at(i + 2);

        // Push fully built triangle onto list_of_triangles vector.
        list_of_triangles.push_back(current_triangle);
      }
    break;

    // Convert list_of_vertices (which represents Quads) into list_of_triangles.

    // NAIVE IMPLEMENTATION of quads. This way of turning quads to triangles is not ideal. There are cases where the triangles pulled from the quad do not represent the original quad. See example below.
    /*
          0-----3
           \     \
            1-----2
    */
    // For instance if we grabbed a,b,c and a,c,d that would be wrong. The program should be made so the user does not have to worry about order of input.

    // UPDATE:
    // The user should expect to start at the first input vertext and draw the next one.
    // The standard convention for describing a shape is from left to right.
    case MGL_QUADS :
      for(unsigned int i = 0; (i + 3) < list_of_vertices.size(); i+=4) {
        triangle current_local_triangle_one;   // first_triangle represents the firt triangle built from the quad.

        // Grab first triangle from quad and store it in current triangle
        current_local_triangle_one.vertex_one   = list_of_vertices.at(i + 0);
        current_local_triangle_one.vertex_two   = list_of_vertices.at(i + 1);
        current_local_triangle_one.vertex_three = list_of_vertices.at(i + 2);

        // Add current_triangle to list of triangles
        list_of_triangles.push_back(current_local_triangle_one);

        triangle current_local_triangle_two;   // first_triangle represents the firt triangle built from the quad.
        // Grab first triangle from quad and store it in current triangle
        current_local_triangle_two.vertex_one   = list_of_vertices.at(i + 2);
        current_local_triangle_two.vertex_two   = list_of_vertices.at(i + 3);
        current_local_triangle_two.vertex_three = list_of_vertices.at(i + 0);

        // Add  second current_triangle to list of triangles
        list_of_triangles.push_back(current_local_triangle_two);

        // // Debugging -- Output Quad Triangle Data
        // cout << "Quad to Triangle Data: " << current_local_triangle_one.vertex_one.pos << endl;
        // cout << "                     : " << current_local_triangle_one.vertex_two.pos << endl;
        // cout << "                     : " << current_local_triangle_one.vertex_three.pos << endl;
        //
        //
        // cout << "Quad to Triangle Data: " << current_local_triangle_two.vertex_one.pos << endl;
        // cout << "                     : " << current_local_triangle_two.vertex_two.pos << endl;
        // cout << "                     : " << current_local_triangle_two.vertex_three.pos << endl;
    }
    break;
  }

  list_of_vertices.clear(); // Clear 'list_of_vertices' after it is converted to 'list_of_triangles'.
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
//  cout << "In mglVertex2" << endl;  //Debugging

  vertex v;

  v.color = current_color;
  v.pos = vec4(x,y,(MGLfloat)0.0,(MGLfloat)1.0);   // Still using a four-dimensional vector. z and w are 0.



  // cout << "\tprojection_matrix():       " << projection_matrix << endl; //Debugging
  // cout << "\tmodelview_matrix():       " << modelview_matrix << endl; //Debugging
  //
  // cout << "\tprojection_matrix() * modelview_matrix():       " << projection_matrix*modelview_matrix << endl; //Debugging
  // cout << "\tv.pos:                     " << v.pos << endl; //Debugging
  //
  // cout << "\tprojection_matrix() * modelview_matrix() * v.pos:       " << projection_matrix*modelview_matrix * v.pos << endl; //Debugging

  mat4 transform_matrix = projection_matrix * modelview_matrix;
  v.pos =  transform_matrix * v.pos;

  //v.pos /= v.pos[3]; // this should occur in mglReadPixels


//  v.pos /= v.pos[3]; // Wrong location


  list_of_vertices.push_back(v);

  // cout << "Final v.pos: " << v.pos << endl; //Debugging
  //cout << "Out mglVertex2" << endl <<endl;  //Debugging
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
//  cout << "In mglVertex3" << endl;  //Debugging
  //cout << z << ", ";
  vertex v;

  v.color = current_color;
  v.pos = vec4(x,y,z,1.0f);

  // cout << "\tprojection_matrix():                      " << projection_matrix << endl; //Debugging
  // cout << "\tmodelview_matrix():                       " << modelview_matrix << endl; //Debugging
  //
  // cout << "\tprojection_matrix() * modelview_matrix(): " << projection_matrix*modelview_matrix << endl; //Debugging
  // cout << "\tv.pos:                     " << v.pos << endl; //Debugging
  //
  // cout << "\tprojection_matrix() * modelview_matrix() * v.pos: " << projection_matrix*modelview_matrix * v.pos << endl; //Debugging

  mat4 transform_matrix = projection_matrix * modelview_matrix;
  v.pos =  transform_matrix * v.pos;

  // v.pos /= v.pos[3]; // Wrong location

  // cout << "Final v.pos: " << v.pos << endl; //Debugging


  list_of_vertices.push_back(v);

  //cout << "Out mglVertex3" << endl <<endl;  //Debugging
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  // Set 'matrix_mode' to the mode passed in by the user.
  matrix_mode = mode;

  // if (mode == MGL_PROJECTION) {
  //   cout << "Current Matrix Mode: projection Mode" << endl;   // Debugging
  // }
  // if (mode == MGL_MODELVIEW) {
  //   cout << "Current Matrix Mode: Model View Mode" << endl;   // Debugging
  // }
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
//  cout << "In mglPushMatrix" << endl; // Debugging
  //cout << "\tget_current_matrix_stack().push_back(get_current_matrix())" << get_current_matrix() <<endl;
  get_current_matrix_stack().push_back(get_current_matrix());
//  cout << "Out mglPushMatrix" << endl << endl;
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
//  cout << "In mglPopMatrix" << endl; // Debugging
  get_current_matrix() = get_current_matrix_stack().back();
  get_current_matrix_stack().pop_back();
//  cout << "Out mglPopMatrix" << endl << endl; // Debugging
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  //cout << "In Load Identity" << endl;
  // Create the identity_matrix
  mat4 identity_matrix;
  identity_matrix.make_zero();     // Start by setting the 'identity_matrix' to 0.

  identity_matrix(0,0)  = (MGLfloat)1;
  identity_matrix(1,1)  = (MGLfloat)1;
  identity_matrix(2,2)  = (MGLfloat)1;
  identity_matrix(3,3)  = (MGLfloat)1;

  get_current_matrix() = identity_matrix;

//  cout << "\tget_current_matrix(): " << get_current_matrix() << endl;

//  cout << "Out Load Identity" << endl << endl;
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
  //cout << "In mglLoadMatrix" << endl;  // Debugging

  for (int i = 0; i < 16; ++i) {
    get_current_matrix().values[i] = *(matrix + i);
  }

//   for (int column = 0; column < 4; column++) {
//     for( int row = 0; row < 4; row++){
//       get_current_matrix().values[row + column] = *(matrix + column + row); // TODO: Need to look at this again
//     }
//   }
// cout << get_current_matrix(); // Debugging
  //cout << "Out mglLoadMatrix" << endl << endl;  // Debugging
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
  //cout << "In mglMultMatrix" << endl; // Debugging
  // // First Entry in matrix
  // get_current_matrix().values[0] = (get_current_matrix().values[0] * *(matrix)) + (get_current_matrix().values[4] * *(matrix + 1)) + (get_current_matrix().values[8] * *(matrix + 2)) + (get_current_matrix().values[12] * *(matrix + 3));
  //
  // // Second Entry in matrix
  // get_current_matrix().values[1] = (get_current_matrix().values[1] * *(matrix)) + (get_current_matrix().values[5] * *(matrix + 1)) + (get_current_matrix().values[9] * *(matrix + 2)) + (get_current_matrix().values[13] * *(matrix + 3));

  mat4 passed_in_matrix;

  for(int i = 0; i < 16; i++) {
    passed_in_matrix.values[i] = *(matrix + i);
  }

  get_current_matrix() = get_current_matrix() * passed_in_matrix;
  // get_current_matrix() =  passed_in_matrix * get_current_matrix();

 //  // Go across columns from a0 towards a12
 //  for(int move_across = 0; move_across < 4; move_across++) {
 //    // Go Down rows from a0 towards a3
 //    for (int i = 0; i < 4; i++) {
 //      // First Entry in matrix
 //      get_current_matrix().values[0 + i]
 //        = (get_current_matrix().values[0 + i] * *(matrix + move_across))
 //        + (get_current_matrix().values[4 + i] * *(matrix + 1 + move_across))
 //        + (get_current_matrix().values[8 + i] * *(matrix + 2 + move_across))
 //        + (get_current_matrix().values[12 + i] * *(matrix + 3 + move_across));
 //    }
 //
 // }


 //cout << "Out mglMultMatrix" << endl << endl;  // Debugging
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
 //cout << "In mglTranslate" << endl; // Debugging

  mat4 translation_matrix;

  translation_matrix.make_zero();

  // Create translation matrix
  translation_matrix(0,0) = 1;
  translation_matrix(0,3) = x;
  translation_matrix(1,1) = 1;
  translation_matrix(1,3) = y;
  translation_matrix(2,2) = 1;
  translation_matrix(2,3) = z;
  translation_matrix(3,3) = 1;


  // cout << "\t x: " << x << ", y: " << y << ", z: " << z << endl;
  // cout << "\tget_current_matrix() = " << get_current_matrix() << endl; // Debugging
  // cout << "\ttranslation_matrix = " << translation_matrix << endl; // Debugging

  get_current_matrix() = get_current_matrix() * translation_matrix;
  // get_current_matrix() = translation_matrix * get_current_matrix();

  // cout << "\tget_current_matrix()*translation_matrix = " << get_current_matrix() << endl; // Debugging

//  cout << "In mglTranslate" << endl << endl; // Debugging
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
  //cout << "In mglRotate" << endl; // Debugging
  mat4 rotation_matrix;

  // convert degrees to radians
  angle = (MGLfloat) (angle*3.14159265358979323846)/180.0;

  MGLfloat magnitude = sqrt(x*x + y*y + z*z);

  if (magnitude != 1) {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  }

  rotation_matrix.make_zero();

  // First row
  rotation_matrix(0,0) = (x*x)*(1-cos(angle)) + cos(angle);    // good
  rotation_matrix(0,1) = (x*y)*(1-cos(angle)) - z*sin(angle);  // good
  rotation_matrix(0,2) = (x*z)*(1-cos(angle)) + y*sin(angle);  // good // first row clear

  // Second row
  rotation_matrix(1,0) = (y*x)*(1-cos(angle)) + z*sin(angle);  // good
  rotation_matrix(1,1) = (y*y)*(1-cos(angle)) + cos(angle);    // good
  rotation_matrix(1,2) = (y*z)*(1-cos(angle)) - x*sin(angle);  // clear

  // Third Row
  rotation_matrix(2,0) = (x*z)*(1-cos(angle)) - y*sin(angle);  // clear
  rotation_matrix(2,1) = (y*z)*(1-cos(angle)) + x*sin(angle);  // clear
  rotation_matrix(2,2) = (z*z)*(1-cos(angle)) + cos(angle);    // clear

  // Fourth row
  rotation_matrix(3,3) = 1;                                    // clear

  //cout << "\t x: " << x << ", y: " << y << ", z: " << z << ", angle: " << angle << ", cos(angle): " << cos(angle) << ", sin(angle): " << sin(angle) << endl;
//  cout << "\tget_current_matrix() = " << get_current_matrix() << endl; // Debugging
  //cout << "\trotation_matrix      = " << rotation_matrix << endl;      // Debugging

  get_current_matrix() = get_current_matrix() * rotation_matrix;
  // get_current_matrix() = rotation_matrix * get_current_matrix();
  //cout << "\tget_current_matrix()*rotation_matrix = " << get_current_matrix() << endl; // Debugging

//  cout << "Out mglRotate" << endl << endl; // Debugging
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  //cout << "In mglScale" << endl; // Debugging
  mat4 scale_matrix;

  scale_matrix.make_zero();

  scale_matrix(0,0) = x;
  scale_matrix(1,1) = y;
  scale_matrix(2,2) = z;
  scale_matrix(3,3) = 1;

  get_current_matrix() = get_current_matrix() * scale_matrix;
  // get_current_matrix() = scale_matrix * get_current_matrix();

//  cout << "In mglScale" << endl << endl; // Debugging
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
//  cout << "In mglFrustum" << endl; // Debugging

  mat4 perspective_matrix;

  perspective_matrix.make_zero();

  perspective_matrix(0,0) = (2.0f*near)/(right-left);      // clear
  perspective_matrix(0,2) = (right+left)/(right-left);     // clear
  perspective_matrix(1,1) = (2.0f*near)/(top-bottom);      // clear
  perspective_matrix(1,2) = (top+bottom)/(top-bottom);     // clear
  perspective_matrix(2,2) = -1.0f*(far+near)/(far-near);   // clear
  perspective_matrix(2,3) = (-2.0f*far*near)/(far-near);   // clear
  perspective_matrix(3,2) = -1.0f;                        // clear

  // cout << "\t left: " << left << ", right: " << right << ", bottom: " << bottom << ", top: " << top << ", near: " << near << ", far: " << far << endl;
  // cout << "\tget_current_matrix() = " << get_current_matrix() << endl; // Debugging
  // cout << "\tperspective_matrix = " << perspective_matrix << endl; // Debugging


  get_current_matrix() = get_current_matrix() * perspective_matrix;
  // get_current_matrix() = perspective_matrix * get_current_matrix();
  // cout << "\tperspective_matrix*get_current_matrix() = " << get_current_matrix() << endl; // Debugging
//  cout << "Out of mglFrustum" << endl << endl; // Debugging
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
//  cout << "In mglOrtho" << endl; // Debugging
  mat4 orthographic_matrix;
  orthographic_matrix.make_zero();

  // Set the diagonal values of the matrix according to the spec on https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/
  orthographic_matrix(0,0) = (MGLfloat) 2/(right - left);
  orthographic_matrix(0,3) = (MGLfloat) -1*((right+left)/(right-left));
  orthographic_matrix(1,1) = (MGLfloat) 2/(top - bottom);
  orthographic_matrix(1,3) = (MGLfloat) -1*((top+bottom)/(top-bottom));
  orthographic_matrix(2,2) = (MGLfloat) -2/(far - near);
  orthographic_matrix(2,3) = (MGLfloat) -1*((far+near)/(far-near));
  orthographic_matrix(3,3) = (MGLfloat) 1;

  get_current_matrix() = get_current_matrix() * orthographic_matrix;
  // get_current_matrix() = orthographic_matrix * get_current_matrix();

  // cout << "\tget_current_matrix() * orthographic_matrix: " << get_current_matrix() << endl; //Debugging
//  cout << "Out of mglOrtho" << endl << endl; // Debugging
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
//  cout << "In mglColor" << endl; // Debugging

  // cout << "\tred: " << red << " green: " << green << " blue: " << blue << endl;

  current_color=vec3(red,green,blue);

//  cout << "Out of mglColor" << endl << endl; // Debugging

}
