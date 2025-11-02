// Transform.cpp: implementation of the Transform class.

// Note: when you construct a matrix using mat4() or mat3(), it will be COLUMN-MAJOR
// Keep this in mind in readfile.cpp and display.cpp
// See FAQ for more details or if you're having problems.

#include "Transform.h"

// Helper rotation function.  Please implement this.  
mat3 Transform::rotate(const float degrees, const vec3& axis) 
{
    // mat3 ret;
    // YOUR CODE FOR HW2 HERE
    // Please implement this.  Likely the same as in HW 1.  
    // return ret;
    float rad = degrees * pi / 180.0;
    mat3 pt1 = mat3(1.0);
    mat3 pt2 = mat3(axis.x * axis.x, axis.x * axis.y, axis.x * axis.z, axis.x * axis.y, axis.y * axis.y, axis.y * axis.z, axis.x * axis.z, axis.y * axis.z, axis.z * axis.z);
    mat3 pt3 = mat3(0, axis.z, -axis.y, -axis.z, 0, axis.x, axis.y, -axis.x, 0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            pt1[i][j] = pt1[i][j] * cos(rad);
            pt2[i][j] = pt2[i][j] * (1 - cos(rad));
            pt3[i][j] = pt3[i][j] * sin(rad);
        }
    }
    //mat3 all = (cos(rad) * pt1) + ((1 - cos(rad)) * pt2) + (sin(rad) * pt3);
    mat3 all = (pt1) + (pt2) + (pt3);
    // You will change this return call
    return all;
}

void Transform::left(float degrees, vec3& eye, vec3& up) 
{
    // YOUR CODE FOR HW2 HERE
    // Likely the same as in HW 1.  
    mat3 rot = rotate(degrees, up);
    eye = rot * eye;
}

void Transform::up(float degrees, vec3& eye, vec3& up) 
{
    // YOUR CODE FOR HW2 HERE 
    // Likely the same as in HW 1.  
    mat3 rot = rotate(degrees, normalize(cross(eye, up)));
    eye = rot * eye;
    up = rot * up;
}

mat4 Transform::lookAt(const vec3 &eye, const vec3 &center, const vec3 &up) 
{
    // mat4 ret;
    // YOUR CODE FOR HW2 HERE
    // Likely the same as in HW 1.  
    // return ret;
    vec3 w = normalize(eye);
    vec3 temp = cross(up, w);
    vec3 u = normalize(temp);
    vec3 v = cross(w, u);
    mat4 r(u.x, v.x, w.x, 0, u.y, v.y, w.y, 0, u.z, v.z, w.z, 0, 0, 0, 0, 1);
    mat4 t(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -eye.x, -eye.y, -eye.z, 1);
    mat4 m = r * t;
    // You will change this return call
    return m;
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar)
{
    float theta = fovy / 2;
    float radians = theta * pi / 180;
    float d = cos(radians) / sin(radians);
    float A = -((zFar + zNear) / (zFar - zNear));
    float B = -((2 * zFar * zNear) / (zFar - zNear));
    mat4 ret = mat4(d / aspect, 0, 0, 0, 0, d, 0, 0, 0, 0, A, -1, 0, 0, B, 0);
    // YOUR CODE FOR HW2 HERE
    // New, to implement the perspective transform as well.  
    // ret = glm::perspective(glm::radians(fovy), aspect, zNear, zFar);
    return ret;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
    mat4 ret = mat4(sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1);
    // YOUR CODE FOR HW2 HERE
    // Implement scaling 
    return ret;
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
    mat4 ret = mat4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, tx, ty, tz, 1);
    // YOUR CODE FOR HW2 HERE
    // Implement translation 
    return ret;
}

// To normalize the up direction and construct a coordinate frame.  
// As discussed in the lecture.  May be relevant to create a properly 
// orthogonal and normalized up. 
// This function is provided as a helper, in case you want to use it. 
// Using this function (in readfile.cpp or display.cpp) is optional.  

vec3 Transform::upvector(const vec3 &up, const vec3 & zvec) 
{
    vec3 x = glm::cross(up,zvec); 
    vec3 y = glm::cross(zvec,x); 
    vec3 ret = glm::normalize(y); 
    return ret; 
}


Transform::Transform()
{

}

Transform::~Transform()
{

}
