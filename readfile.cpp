/*****************************************************************************/
/* This is the program skeleton for homework 2 in CSE167 by Ravi Ramamoorthi */
/* Extends HW 1 to deal with shading, more transforms and multiple objects   */
/*****************************************************************************/

/*****************************************************************************/
// This file is readfile.cpp.  It includes helper functions for matrix 
// transformations for a stack (matransform) and to rightmultiply the 
// top of a stack.  These functions are given to aid in setting up the 
// transformations properly, and to use glm functions in the right way.  
// Their use is optional in your program.  
  

// The functions readvals and readfile do basic parsing.  You can of course 
// rewrite the parser as you wish, but we think this basic form might be 
// useful to you.  It is a very simple parser.

// Please fill in parts that say YOUR CODE FOR HW 2 HERE. 
// Read the other parts to get a context of what is going on. 
  
/*****************************************************************************/

// Basic includes to get this file to work.  
#include <iostream>
#include <string>
#include <fstream>
//#include <sstream>
#include <deque>
//#include <stack>

#include "Transform.h" 

using namespace std;
#include "variables.h" 
#include "readfile.h"

// You may not need to use the following two functions, but it is provided
// here for convenience
void rightmultiply(const mat4 & M, stack<mat4> &transfstack) 
{
    mat4 &T = transfstack.top(); 
    T = T * M; 
}

bool readvals(stringstream &s, const int numvals, float* values) 
{
    for (int i = 0; i < numvals; i++) {
        s >> values[i]; 
        if (s.fail()) {
            cout << "Failed reading value " << i << " will skip\n"; 
            return false;
        }
    }
    return true; 
}

// The function below applies the appropriate transform to a 4-vector
void readfile(const char* filename) 
{
    string str, cmd; 
    ifstream in;
    in.open(filename); 
    if (in.is_open()) {

        // I need to implement a matrix stack to store transforms.  
        // This is done using standard STL Templates 
        stack <mat4> transfstack; 
        transfstack.push(mat4(1.0));  // identity

        getline (in, str); 
        while (in) {
            if ((str.find_first_not_of(" \t\r\n") != string::npos) && (str[0] != '#')) {
                // Ruled out comment and blank lines 

                stringstream s(str);
                s >> cmd; 
                int i; 
                float values[10]; // Position and color for light, colors for others
                                    // Up to 10 params for cameras.  
                bool validinput; // Validity of input 

                // Material Commands 
                // Ambient, diffuse, specular, shininess properties for each object.
                // Filling this in is pretty straightforward, so I've left it in 
                // the skeleton, also as a hint of how to do the more complex ones.
                // Note that no transforms/stacks are applied to the colors. 

                if (cmd == "maxverts") {
                    validinput = readvals(s, 1, values);
                    if (validinput) {
                        maxverts = values[0];
                    }
                } else if (cmd == "maxvertnorms") {
                    validinput = readvals(s, 1, values);
                    if (validinput) {
                        maxvertnorms = values[0];
                    }
                } else if (cmd == "directional") {
                    validinput = readvals(s, 6, values); 
                    if (validinput) {
                        for (i = 0; i < 6; i++) {
                            directional[i] = values[i]; 
                            dirExits = 1;
                        }
                    }
                } else if (cmd == "point") {
                    validinput = readvals(s, 6, values); 
                    if (validinput) {
                        for (i = 0; i < 6; i++) {
                            point[i] = values[i]; 
                            pntExits = 1;
                        }
                    }
                } else if (cmd == "attenuation") {
                    validinput = readvals(s, 3, values); 
                    if (validinput) {
                        attenExists = 1;
                        for (i = 0; i < 3; i++) {
                            attenuation[i] = values[i]; 
                        }
                    }
                } else if (cmd == "ambient") {
                    validinput = readvals(s, 3, values); // colors 
                    if (validinput) {
                        for (i = 0; i < 3; i++) {
                            ambient[i] = values[i]; 
                        }
                    }
                } else if (cmd == "diffuse") {
                    validinput = readvals(s, 3, values); // colors 
                    if (validinput) {
                        for (i = 0; i < 3; i++) {
                            diffuse[i] = values[i]; 
                        }
                    }
                } else if (cmd == "specular") {
                    validinput = readvals(s, 3, values); // colors 
                    if (validinput) {
                        for (i = 0; i < 3; i++) {
                            specular[i] = values[i]; 
                        }
                    }
                } else if (cmd == "shininess") {
                    validinput = readvals(s, 1, values); // colors 
                    if (validinput) {
                        shininess = values[0];
                    }
                } else if (cmd == "emission") {
                    validinput = readvals(s, 3, values); // colors 
                    if (validinput) {
                        for (i = 0; i < 3; i++) {
                            emission[i] = values[i]; 
                        }
                    }
                } else if (cmd == "size") {
                    validinput = readvals(s,2,values); 
                    if (validinput) { 
                        w = (int) values[0]; h = (int) values[1]; 
                        width = values[0];
                        height = values[1];
                    } 
                } else if (cmd == "maxdepth") {
                    validinput = readvals(s,1,values); 
                    if (validinput) { 
                        depth = (int) values[0]; 
                    } 
                } else if (cmd == "camera") {
                    validinput = readvals(s,10,values); // 10 values eye cen up fov
                    if (validinput) {

                        // YOUR CODE FOR HW 2 HERE
                        // Use all of values[0...9]
                        // You may need to use the upvector fn in Transform.cpp
                        // to set up correctly. 
                        // Set eyeinit upinit center fovy in variables.h 
                        eyeinit = vec3(values[0], values[1], values[2]);
                        center = vec3(values[3], values[4], values[5]);
                        upinit = normalize(vec3(values[6], values[7], values[8]));
                        fovy = values[9];
                        yfov = values[9] * pi / 180.0f;
                        xfov = 2.0f * atan(tan(yfov / 2.0f) * ((float)width / (float)height));
                        //fovx = 2 * atan(tan(((fovy * (pi / 180)) / 2) * (w / h)));
                    }
                }

                // I've left the code for loading objects in the skeleton, so 
                // you can get a sense of how this works.  
                // Also look at demo.txt to get a sense of why things are done this way.
                else if (cmd == "vertex") {
                    if (numverts == maxverts) {
                        cerr << "Reached Maximum Number of vertices " << numverts << "Will ignore further vertices\n";
                    } else {
                        validinput = readvals(s, 3, values);
                        if (validinput) {
                            vertex * vert = &(vertexpile[numverts]);
                            vert->v.x = values[0];
                            vert->v.y = values[1];
                            vert->v.z = values[2];
                            numverts++;
                        }
                    }
                } else if (cmd == "vertexnormal") {
                    if (numverts == maxverts) {
                        cerr << "Reached Maximum Number of vertices " << numverts << "Will ignore further vertices\n";
                    } else {
                        validinput = readvals(s, 6, values);
                        if (validinput) {
                            vertexnormal * vert = &(vertexnormalpile[numverts]);
                            vert->x = values[0];
                            vert->y = values[1];
                            vert->z = values[2];
                            vert->nx = values[3];
                            vert->ny = values[4];
                            vert->nz = values[5];
                            numvertnorms++;
                        }
                    }
                } else if (cmd == "tri") {
                    validinput = readvals(s, 3, values);
                    if (validinput) {
                        triangle * tri = &(triangles[numtriangles]);
                        for (int i = 0; i < 3; i++) {
                            tri->ambient[i] = ambient[i];
                            tri->diffuse[i] = diffuse[i];
                            tri->specular[i] = specular[i];
                            tri->emission[i] = emission[i];
                        }
                        tri->shininess = shininess;
                        tri->v1 = vertexpile[(int)values[0]].v;
                        tri->v2 = vertexpile[(int)values[1]].v;
                        tri->v3 = vertexpile[(int)values[2]].v;
                        tri->trans = transfstack.top();
                        numtriangles++;
                    }
                } else if (cmd == "trinormal") {
                    validinput = readvals(s, 3, values);
                    if (validinput) {
                        trianglenormal * tri = &(trianglenormals[numtrianglenormals]);
                        for (int i = 0; i < 3; i++) {
                            tri->ambient[i] = ambient[i];
                            tri->diffuse[i] = diffuse[i];
                            tri->specular[i] = specular[i];
                            tri->emission[i] = emission[i];
                        }
                        tri->shininess = shininess;
                        tri->v1 = vertexnormalpile[(int)values[0]];
                        tri->v2 = vertexnormalpile[(int)values[1]];
                        tri->v3 = vertexnormalpile[(int)values[2]];
                        numtrianglenormals++;
                    }
                } else if (cmd == "sphere") {
                    validinput = readvals(s, 4, values);
                    if (validinput) {
                        sph * sp = &(spheres[numspheres]);
                        for (int i = 0; i < 3; i++) {
                            sp->ambient[i] = ambient[i];
                            sp->diffuse[i] = diffuse[i];
                            sp->specular[i] = specular[i];
                            sp->emission[i] = emission[i];
                        }
                        sp->shininess = shininess;
                        sp->v.x = values[0];
                        sp->v.y = values[1];
                        sp->v.z = values[2];
                        sp->trans = transfstack.top();
                        sp->radius = values[3];
                        numspheres++;
                    }
                }
                /*
                else if (cmd == "sphere" || cmd == "cube" || cmd == "teapot") {
                    if (numobjects == maxobjects) { // No more objects 
                        cerr << "Reached Maximum Number of Objects " << numobjects << " Will ignore further objects\n";
                    } else {
                        validinput = readvals(s, 1, values); 
                        if (validinput) {
                            object * obj = &(objects[numobjects]); 
                            obj->size = values[0]; 

                            // Set the object's light properties
                            for (i = 0; i < 4; i++) {
                                (obj->ambient)[i] = ambient[i]; 
                                (obj->diffuse)[i] = diffuse[i]; 
                                (obj->specular)[i] = specular[i]; 
                                (obj->emission)[i] = emission[i];
                            }
                            obj->shininess = shininess; 

                            // Set the object's transform
                            obj->transform = transfstack.top(); 

                            // Set the object's type
                            if (cmd == "sphere") {
                                obj->type = sphere; 
                            } else if (cmd == "cube") {
                                obj->type = cube; 
                            } else if (cmd == "teapot") {
                                obj->type = teapot; 
                            }
                        }
                        ++numobjects; 
                    }
                }
                */
                else if (cmd == "translate") {
                    validinput = readvals(s,3,values); 
                    if (validinput) {

                        // YOUR CODE FOR HW 2 HERE.  
                        // Think about how the transformation stack is affected
                        // You might want to use helper functions on top of file. 
                        // Also keep in mind what order your matrix is!
                        mat4 M = Transform::translate(values[0], values[1], values[2]);

         

                        rightmultiply(M, transfstack);
                    }
                }
                else if (cmd == "scale") {
                    validinput = readvals(s,3,values); 
                    if (validinput) {

                        // YOUR CODE FOR HW 2 HERE.  
                        // Think about how the transformation stack is affected
                        // You might want to use helper functions on top of file.  
                        // Also keep in mind what order your matrix is!
                        mat4 M = Transform::scale(values[0], values[1], values[2]);

                      

                        rightmultiply(M, transfstack);
              
                    }
                }
                else if (cmd == "rotate") {
                    validinput = readvals(s,4,values); 
                    if (validinput) {

                        // YOUR CODE FOR HW 2 HERE. 
                        // values[0..2] are the axis, values[3] is the angle.  
                        // You may want to normalize the axis (or in Transform::rotate)
                        // See how the stack is affected, as above.  
                        // Note that rotate returns a mat3. 
                        // Also keep in mind what order your matrix is!

                        vec3 axis = vec3(values[0], values[1], values[2]);

                        mat3 matrix3 = Transform::rotate(values[3], axis);

                        mat4 M = mat4(matrix3[0][0], matrix3[0][1], matrix3[0][2], 0, matrix3[1][0], matrix3[1][1], matrix3[1][2], 0, matrix3[2][0], matrix3[2][1], matrix3[2][2], 0, 0, 0, 0, 1);

                  

                        rightmultiply(M, transfstack);

                    }
                }

                // I include the basic push/pop code for matrix stacks
                else if (cmd == "pushTransform") {
                    transfstack.push(transfstack.top()); 
                } else if (cmd == "popTransform") {
                    if (transfstack.size() <= 1) {
                        cerr << "Stack has no elements.  Cannot Pop\n"; 
                    } else {
                        transfstack.pop(); 
                    }
                }

                else {
                    cerr << "Unknown Command: " << cmd << " Skipping \n"; 
                }
            }
            /*
            for (int i = 0; i < 40; i++) {
                cout << "Value of i is " << i << " Value of ambient is " << ambient[i] << "\n";
            }
            */
            getline (in, str); 
        }

        // Set up initial position for eye, up and amount
        // As well as booleans 

        eye = eyeinit; 
        up = upinit; 
        amount = 5;
        sx = sy = 1.0;  // keyboard controlled scales in x and y 
        tx = ty = 0.0;  // keyboard controllled translation in x and y  

    } else {
        cerr << "Unable to Open Input Data File " << filename << "\n"; 
        throw 2; 
    }
}
