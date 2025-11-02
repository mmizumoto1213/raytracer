/******************************************************************************/
/* This is the program skeleton for homework 2 in CSE167 by Ravi Ramamoorthi  */
/* Extends HW 1 to deal with shading, more transforms and multiple objects    */
/******************************************************************************/

// You shouldn't have to edit this file at all!

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <deque>
#include <stack>
#include <omp.h>
#include "Transform.h"
#include <C:\Users\mmizu\OneDrive\Desktop\hw2-windows\FreeImage.h>

//#include <C:\Users\mmizu\OneDrive\Desktop\School\CSE167\hw4-windows\hw2-windows\hw2-windows\FreeImage.lib>
//#include <C:\Users\mmizu\OneDrive\Desktop\School\CSE167\hw4-windows\hw2-windows\hw2-windows\FreeImage.dll>

// Main variables in the program.  
#define MAINPROGRAM 
#include "variables.h" 
#include "readfile.h" // prototypes for readfile.cpp  

class ray {
  public:
  ray() {}

  ray(const vec3& origin, const vec3& direction) {
    ori = origin;
    dir = direction;
  }

  vec3 origin() const  { return ori; }
  vec3 direction() const { return dir; }

  private:
    vec3 ori;
    vec3 dir;
};

class intersection {
  public:
  intersection() {
    n = 0;
    mind = INFINITY;
  }

  intersection(const vec3& position, const vec3& normal, const float mindist, const int val) {
    n = 1;
    pos = position;
    norm = normal;
    mind = mindist;
    type = val;
  }

  int getN() const { return n; }
  vec3 getPosition() const  { return pos; }
  vec3 getNormal() const  { return norm; }
  float getDistance() const  { return mind; }
  int getShape() const { return type; }
  triangle getTriangle() const { return tri; }
  sph getSphere() const { return sphere; }
  void setTriangle(triangle t) {
    tri = t;
  }
  void setSphere(sph s) {
    sphere = s;
  }

  private:
    int n;  // 1 if intersection exists 0 if doesnt
    vec3 pos;
    vec3 norm;
    float mind;
    int type; // 1 means triangle 2 means sphere
    triangle tri;
    sph sphere;
};

class object {
  public:
  object() {}

  object(int t) {
    type = t;
  }

  int getType() { return type; }
  triangle getTriangle() { return tri; }
  sph getSphere() { return sphere; }
  void setTriangle(triangle tr) {
    tri = tr;
  }
  void setSphere(sph sp) {
    sphere = sp;
  }

  private:
    int type; // 1 means triangle 2 means sphere
    triangle tri;
    sph sphere;
};

class bbox {
  public:

  bbox() {
    type = 0;
    leftexists = 0;
    rightexists = 0;
  }

  bbox(int t) {
    type = t;
    leftexists = 0;
    rightexists = 0;
  }

  bbox(const vec3& vbig, const vec3& vsmall) {
    type = 0;
    v1 = vbig;
    v2 = vsmall;
  }

  vec3 getVBig() { return v1; }
  vec3 getVSmall() { return v2; }
  bbox getLeft() { return *left; }
  bbox getRight() { return *right; }
  int getLeftExists() { return leftexists; }
  int getRightExists() { return rightexists; }
  int getType() { return type; }
  sph getSphere() { return sphere; }
  triangle getTriangle() { return tri; }
  void setVBig(vec3 big) {
    v1 = big;
  }
  void setVSmall(vec3 small) {
    v2 = small;
  }
  void setLeft(bbox *l) {
    left = l;
    leftexists = 1;
  } 
  void setRight(bbox *r) {
    right = r;
    rightexists = 1;
  }
  void setTriangle(triangle t) {
    tri = t;
  }
  void setSphere(sph s) {
    sphere = s;
  }

  private:
    int type; // 1 means triangle 2 means sphere 3 means bbox
    vec3 v1;
    vec3 v2;
    struct bbox *left;
    struct bbox *right;
    int leftexists;
    int rightexists;
    triangle tri;
    sph sphere;
};

bbox bboxcombine (bbox b1, bbox b2) {
  float maxx = max(b1.getVBig().x, b2.getVBig().x);
  float maxy = max(b1.getVBig().y, b2.getVBig().y);
  float maxz = max(b1.getVBig().z, b2.getVBig().z);
  float minx = min(b1.getVSmall().x, b2.getVSmall().x);
  float miny = min(b1.getVSmall().y, b2.getVSmall().y);
  float minz = min(b1.getVSmall().z, b2.getVSmall().z);
  bbox toRet;
  toRet.setVBig(vec3(maxx, maxy, maxz));
  toRet.setVSmall(vec3(minx, miny, minz));
  return toRet;
}

bbox bboxbound (object a) {
  bbox b;
  if (a.getType() == 1) {
    float maxx = a.getTriangle().v1.x;
    float minx = a.getTriangle().v1.x;
    float maxy = a.getTriangle().v1.y;
    float miny = a.getTriangle().v1.y;
    float maxz = a.getTriangle().v1.z;
    float minz = a.getTriangle().v1.z;
    if (a.getTriangle().v2.x > maxx) {
      float maxx = a.getTriangle().v2.x;
    } else if (a.getTriangle().v2.y > maxy) {
      float maxy = a.getTriangle().v2.y;
    } else if (a.getTriangle().v2.z > maxz) {
      float maxz = a.getTriangle().v2.z;
    }
    if (a.getTriangle().v3.x > maxx) {
      float maxx = a.getTriangle().v3.x;
    } else if (a.getTriangle().v3.y > maxy) {
      float maxy = a.getTriangle().v3.y;
    } else if (a.getTriangle().v3.z > maxz) {
      float maxz = a.getTriangle().v3.z;
    }
    if (a.getTriangle().v2.x < minx) {
      float minx = a.getTriangle().v2.x;
    } else if (a.getTriangle().v2.y < miny) {
      float miny = a.getTriangle().v2.y;
    } else if (a.getTriangle().v2.z < minz) {
      float minz = a.getTriangle().v2.z;
    }
    if (a.getTriangle().v3.x < minx) {
      float minx = a.getTriangle().v3.x;
    } else if (a.getTriangle().v3.y < miny) {
      float miny = a.getTriangle().v3.y;
    } else if (a.getTriangle().v3.z < minz) {
      float minz = a.getTriangle().v3.z;
    }
    b.setVBig(vec3(maxx, maxy, maxz));
    b.setVSmall(vec3(minx, miny, minz));
    //b.o = a[0];
  } else if (a.getType() == 2) {
    float maxx = a.getSphere().v.x + a.getSphere().radius;
    float maxy = a.getSphere().v.y + a.getSphere().radius;
    float maxz = a.getSphere().v.z + a.getSphere().radius;
    float minx = a.getSphere().v.x - a.getSphere().radius;
    float miny = a.getSphere().v.y - a.getSphere().radius;
    float minz = a.getSphere().v.z - a.getSphere().radius;
    b.setVBig(vec3(maxx, maxy, maxz));
    b.setVSmall(vec3(minx, miny, minz));
    //b.o = a[0];
  } 
  return b;
}

bbox bboxcreate (object a[pilesize], int size, int axis) {
  bbox b;
  if (size == 1) {
    bboxbound (a[0]);
    if (a[0].getType() == 1) {
      bbox left = bbox(1);
      left.setTriangle(a[0].getTriangle());
      b.setLeft(&left);
    } else {
      bbox left = bbox(2);
      left.setSphere(a[0].getSphere());
      b.setLeft(&left);
    }
    return b;
  } else if (size == 2) {
    bbox left = bboxbound(a[0]);
    bbox right = bboxbound(a[1]);
    b = bboxcombine(left, right);
    b.setLeft(&left);
    b.setRight(&right);
    return b;
  } else {
    // axis = 0 means x axis, 1 means y axis, 2 means z axis
    if (axis == 0) {
      float maxx;
      float minx;
      if (a[0].getType() == 1) {
        maxx = a[0].getTriangle().v1.x;
        minx = maxx;
      } else if (a[0].getType() == 2) {
        maxx = a[0].getSphere().v.x;
        minx = maxx;
      } 
      for (int i = 0; i < size; i++) {
        if (a[i].getType() == 1) {
          maxx = max(a[i].getTriangle().v1.x, maxx);
          maxx = max(a[i].getTriangle().v2.x, maxx);
          maxx = max(a[i].getTriangle().v3.x, maxx);
          minx = min(a[i].getTriangle().v1.x, minx);
          minx = min(a[i].getTriangle().v2.x, minx);
          minx = min(a[i].getTriangle().v3.x, minx);
        } else if (a[i].getType() == 2) {
          maxx = max(a[i].getSphere().v.x, maxx);
          minx = min(a[i].getSphere().v.x, minx);
        }
      } 
      float middlex = (maxx + minx) / 2.0f;
      int leftcount = 0;
      int rightcount = 0;
      object lft[pilesize];
      object rgt[pilesize];
      for (int i = 0; i < size; i++) {
        if (a[i].getType() == 1) {
          if (a[i].getTriangle().v1.x < middlex) {
            lft[leftcount] = a[i];
            leftcount++;
          } else {
            rgt[rightcount] = a[i];
            rightcount++;
          }
        } else if (a[i].getType() == 2) {
          if (a[i].getSphere().v.x < middlex) {
            lft[leftcount] = a[i];
            leftcount++;
          } else {
            rgt[rightcount] = a[i];
            rightcount++;
          }
        } 
      }
      bbox left = bboxcreate (lft, leftcount, (axis + 1) % 3);
      bbox right = bboxcreate (rgt, rightcount, (axis + 1) % 3);
      b = bboxcombine(left, right);
      b.setLeft(&left);
      b.setRight(&right);
    } else if (axis == 1) {
      float maxy;
      float miny;
      if (a[0].getType() == 1) {
        maxy = a[0].getTriangle().v1.y;
        miny = maxy;
      } else if (a[0].getType() == 2) {
        maxy = a[0].getSphere().v.y;
        miny = maxy;
      } 
      for (int i = 0; i < size; i++) {
        if (a[i].getType() == 1) {
          maxy = max(a[i].getTriangle().v1.y, maxy);
          maxy = max(a[i].getTriangle().v2.y, maxy);
          maxy = max(a[i].getTriangle().v3.y, maxy);
          miny = min(a[i].getTriangle().v1.y, miny);
          miny = min(a[i].getTriangle().v2.y, miny);
          miny = min(a[i].getTriangle().v3.y, miny);
        } else if (a[i].getType() == 2) {
          maxy = max(a[i].getSphere().v.y, maxy);
          miny = min(a[i].getSphere().v.y, miny);
        } 
      } 
      float middley = (maxy + miny) / 2.0f;
      int leftcount = 0;
      int rightcount = 0;
      object lft[pilesize];
      object rgt[pilesize];
      for (int i = 0; i < size; i++) {
        if (a[i].getType() == 1) {
          if (a[i].getTriangle().v1.y < middley) {
            lft[leftcount] = a[i];
            leftcount++;
          } else {
            rgt[rightcount] = a[i];
            rightcount++;
          }
        } else if (a[i].getType() == 2) {
          if (a[i].getSphere().v.y < middley) {
            lft[leftcount] = a[i];
            leftcount++;
          } else {
            rgt[rightcount] = a[i];
            rightcount++;
          }
        } 
      }
      bbox left = bboxcreate (lft, leftcount, (axis + 1) % 3);
      bbox right = bboxcreate (rgt, rightcount, (axis + 1) % 3);
      b = bboxcombine(left, right);
      b.setLeft(&left);
      b.setRight(&right);
    } else if (axis == 2) {
      float maxz;
      float minz;
      if (a[0].getType() == 1) {
        maxz = a[0].getTriangle().v1.z;
        minz = maxz;
      } else if (a[0].getType() == 2) {
        maxz = a[0].getSphere().v.z;
        minz = maxz;
      } 
      for (int i = 0; i < size; i++) {
        if (a[i].getType() == 1) {
          maxz = max(a[i].getTriangle().v1.z, maxz);
          maxz = max(a[i].getTriangle().v2.z, maxz);
          maxz = max(a[i].getTriangle().v3.z, maxz);
          minz = min(a[i].getTriangle().v1.z, minz);
          minz = min(a[i].getTriangle().v2.z, minz);
          minz = min(a[i].getTriangle().v3.z, minz);
        } else if (a[i].getType() == 2) {
          maxz = max(a[i].getSphere().v.z, maxz);
          minz = min(a[i].getSphere().v.z, minz);
        } 
      } 
      float middlex = (maxz + minz) / 2.0f;
      int leftcount = 0;
      int rightcount = 0;
      object lft[pilesize];
      object rgt[pilesize];
      for (int i = 0; i < size; i++) {
        if (a[i].getType() == 1) {
          if (a[i].getTriangle().v1.z < middlex) {
            lft[leftcount] = a[i];
            leftcount++;
          } else {
            rgt[rightcount] = a[i];
            rightcount++;
          }
        } else if (a[i].getType() == 2) {
          if (a[i].getSphere().v.z < middlex) {
            lft[leftcount] = a[i];
            leftcount++;
          } else {
            rgt[rightcount] = a[i];
            rightcount++;
          }
        } 
      }
      bbox left = bboxcreate (lft, leftcount, (axis + 1) % 3);
      bbox right = bboxcreate (rgt, rightcount, (axis + 1) % 3);
      b = bboxcombine(left, right);
      b.setLeft(&left);
      b.setRight(&right);
    } 
    return b;
  }
}





intersection IntersectBBOX (const ray r, bbox b) {
  intersection intersec;
  if (b.getType() == 0) { // BBOX
    float txmin = (b.getVSmall().x - r.origin().x) / (r.direction().x);
    float txmax = (b.getVBig().x - r.origin().x) / (r.direction().x);
    float tymin = (b.getVSmall().y - r.origin().y) / (r.direction().y);
    float tymax = (b.getVBig().y - r.origin().y) / (r.direction().y);
    float tzmin = (b.getVSmall().z - r.origin().z) / (r.direction().z);
    float tzmax = (b.getVBig().z - r.origin().z) / (r.direction().z);
    if (txmin > tymax || txmin > tzmax || tymin > txmax || tymin > tzmax || tzmin > txmax || tzmin > tymax) {
      // no intersection
    } else {
      if (b.getLeftExists()) {
        intersec = IntersectBBOX (r, b.getLeft());
      }
      if (b.getRightExists()) {
        intersec = IntersectBBOX (r, b.getRight());
      }
    }
  } else if (b.getType() == 2) { // SPHERE
    sph sphere = b.getSphere();
    vec3 cent = sphere.v;
    mat4 inv = inverse(sphere.trans);
    vec4 centtranstemp = inv * vec4(cent, 1.0f);
    vec4 rayoritemp = inv * vec4(r.origin(), 1.0f);
    vec4 raydirtemp = inv * vec4(r.direction(), 0.0f);
    vec3 centtrans = vec3(centtranstemp.x, centtranstemp.y, centtranstemp.z);
    vec3 rayori = vec3(rayoritemp.x, rayoritemp.y, rayoritemp.z);
    vec3 raydir = vec3(raydirtemp.x, raydirtemp.y, raydirtemp.z);
    ray raytrans = ray (rayori, raydir);
    // float a = dot(r.direction() , r.direction());
    // float b = 2 * dot(r.direction(), r.origin() - cent);
    // float c = dot(r.origin() - cent, r.origin() - cent) - (sphere.radius * sphere.radius);
    float a = dot(raydir , raydir);
    float b = 2 * dot(raydir, rayori - cent);
    float c = dot(rayori - cent, rayori - cent) - (sphere.radius * sphere.radius);
    float discriminant = (b * b) - (4 * a * c);
    //cout << "raydirection: " << " " << r.direction().x << " " << r.direction().y << " " <<r.direction().z << '\n';
    //cout << "a = " << a << " b = " << b << " c = " << c << endl;
    //cout << discriminant << endl;
    if (discriminant >= 0) {
      //cout << "potential" << endl;
      float t1 = (-b + sqrt(discriminant)) / (2 * a);
      float t2 = (-b - sqrt(discriminant)) / (2 * a);
      // Case both pos take the smaller of the two;
      if (t1 > 0 && t2 > 0) {
        //cout << "potential" << endl;
        float t = min (t1, t2);
        if (t < intersec.getDistance()) {
          // mindist = t;
          vec3 point = rayori + (raydir * t);
          vec4 pointtemp = vec4(point, 1.0f);
          vec4 pointtemp2 = sphere.trans * pointtemp;
          vec3 pointposition = vec3(pointtemp2.x, pointtemp2.y, pointtemp2.z);
          // vec3 ntransnorm = r.origin() + (r.direction() * mindist);
          vec4 normtemp = inverseTranspose(sphere.trans) * vec4((point), 1.0f);
          vec3 norm = normalize (vec3(normtemp.x, normtemp.y, normtemp.z)) * -1.0f;
          // vec3 normpoint = r.origin() + (r.direction() * mindist) - cent;
          // vec4 normpointtrans = inv * vec4(normpoint, 1.0f);
          // vec3 norm = normalize (vec3(normpointtrans.x, normpointtrans.y, normpointtrans.z)) * -1.0f;
          intersec = intersection(pointposition, norm, t, 2);
          intersec.setSphere(sphere);
        }
      // Case that one is pos
      } else if (t1 > 0 || t2 > 0) {
        //cout << "potential" << endl;
        float t = max (t1, t2);
        if (t < intersec.getDistance()) {
          // mindist = t;
          vec3 point = rayori + (raydir * t);
          vec4 pointtemp = vec4(point, 1.0f);
          vec4 pointtemp2 = sphere.trans * pointtemp;
          vec3 pointposition = vec3(pointtemp2.x, pointtemp2.y, pointtemp2.z);
          // vec3 ntransnorm = r.origin() + (r.direction() * mindist);
          vec4 normtemp = inverseTranspose(sphere.trans) * vec4((point), 1.0f);
          vec3 norm = normalize (vec3(normtemp.x, normtemp.y, normtemp.z)) * -1.0f;
          // vec3 normpoint = r.origin() + (r.direction() * mindist) - cent;
          // vec4 normpointtrans = inv * vec4(normpoint, 1.0f);
          // vec3 norm = normalize (vec3(normpointtrans.x, normpointtrans.y, normpointtrans.z)) * -1.0f;
          intersec = intersection(pointposition, norm, t, 2);
          intersec.setSphere(sphere);
        }
      }
      //cout << "t1: " << t1 << " t2: " << t2 << endl;
    }
  } else if (b.getType() == 1) { // TRIANGLE
    triangle tri = b.getTriangle();
    vec4 tempa = tri.trans * vec4(tri.v1, 1.0f);
    vec4 tempb = tri.trans * vec4(tri.v2, 1.0f);
    vec4 tempc = tri.trans * vec4(tri.v3, 1.0f);
    // vec3 a = tri.v1;
    // vec3 b = tri.v2;
    // vec3 c = tri.v3;
    vec3 a = vec3(tempa.x, tempa.y, tempa.z);
    vec3 b = vec3(tempb.x, tempb.y, tempb.z);
    vec3 c = vec3(tempc.x, tempc.y, tempc.z);
    vec3 norm = normalize(cross(c-a, b-a));
    float denomcheck = dot(r.direction(), norm);
    if (denomcheck != 0) {
      float t = ((dot(a, norm) - dot(r.origin(), norm)) / dot(r.direction(), norm));
      vec3 pointp = r.origin() + (r.direction() * t);
      vec3 pminusa = pointp - a;
      vec3 bminusa = b - a;
      vec3 cminusa = c - a;
      float gammanumerator = (pminusa.y * bminusa.x) - (pminusa.x * bminusa.y);
      float gammadenominator = (bminusa.y * cminusa.x) - (cminusa.y * bminusa.x);
      float gamma = -1.0f * (gammanumerator / gammadenominator);
      float betapt1 = (pminusa.x / bminusa.x);
      float betapt2 = (gamma * cminusa.x) / bminusa.x;
      float beta = betapt1 - betapt2;
      float alpha = (pointp.x - (beta * b.x) - (gamma * c.x)) / a.x;
      float total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < intersec.getDistance()) || (t > 0 && intersec.getDistance() == -1)) {
          // mindist = t;

          intersec = intersection((r.origin() + (r.direction() * t)), norm, t, 1);
          intersec.setTriangle(tri);
          //cout << "tri found" << endl;
        }
      } 
      gammanumerator = (pminusa.z * bminusa.y) - (pminusa.y * bminusa.z);
      gammadenominator = (bminusa.z * cminusa.y) - (cminusa.z * bminusa.y);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.y / bminusa.y);
      betapt2 = (gamma * cminusa.y) / bminusa.y;
      beta = betapt1 - betapt2;
      alpha = (pointp.y - (beta * b.y) - (gamma * c.y)) / a.y;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < intersec.getDistance()) || (t > 0 && intersec.getDistance() == -1)) {
          // mindist = t;

          intersec = intersection((r.origin() + (r.direction() * t)), norm, t, 1);
          intersec.setTriangle(tri);
          //cout << "tri found" << endl;
        }
      }
      gammanumerator = (pminusa.x * bminusa.z) - (pminusa.z * bminusa.x);
      gammadenominator = (bminusa.x * cminusa.z) - (cminusa.x * bminusa.z);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.z / bminusa.z);
      betapt2 = (gamma * cminusa.z) / bminusa.z;
      beta = betapt1 - betapt2;
      alpha = (pointp.z - (beta * b.z) - (gamma * c.z)) / a.z;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < intersec.getDistance()) || (t > 0 && intersec.getDistance() == -1)) {
          // mindist = t;

          intersec = intersection((r.origin() + (r.direction() * t)), norm, t, 1);
          intersec.setTriangle(tri);
          //cout << "tri found" << endl;
        }
      }
      gammanumerator = (pminusa.x * bminusa.y) - (pminusa.y * bminusa.x);
      gammadenominator = (bminusa.x * cminusa.y) - (cminusa.x * bminusa.y);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.y / bminusa.y);
      betapt2 = (gamma * cminusa.y) / bminusa.y;
      beta = betapt1 - betapt2;
      alpha = (pointp.y - (beta * b.y) - (gamma * c.y)) / a.y;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < intersec.getDistance()) || (t > 0 && intersec.getDistance() == -1)) {
          // mindist = t;

          intersec = intersection((r.origin() + (r.direction() * t)), norm, t, 1);
          intersec.setTriangle(tri);
          //cout << "tri found" << endl;
        }
      } 
      gammanumerator = (pminusa.y * bminusa.z) - (pminusa.z * bminusa.y);
      gammadenominator = (bminusa.y * cminusa.z) - (cminusa.y * bminusa.z);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.z / bminusa.z);
      betapt2 = (gamma * cminusa.z) / bminusa.z;
      beta = betapt1 - betapt2;
      alpha = (pointp.z - (beta * b.z) - (gamma * c.z)) / a.z;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < intersec.getDistance()) || (t > 0 && intersec.getDistance() == -1)) {
          // mindist = t;

          intersec = intersection((r.origin() + (r.direction() * t)), norm, t, 1);
          intersec.setTriangle(tri);
          //cout << "tri found" << endl;
        }
      }
      gammanumerator = (pminusa.z * bminusa.x) - (pminusa.x * bminusa.z);
      gammadenominator = (bminusa.z * cminusa.x) - (cminusa.z * bminusa.x);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.x / bminusa.x);
      betapt2 = (gamma * cminusa.x) / bminusa.x;
      beta = betapt1 - betapt2;
      alpha = (pointp.x - (beta * b.x) - (gamma * c.x)) / a.x;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < intersec.getDistance()) || (t > 0 && intersec.getDistance() == -1)) {
          // mindist = t;

          intersec = intersection((r.origin() + (r.direction() * t)), norm, t, 1);
          intersec.setTriangle(tri);
          //cout << "tri found" << endl;
        }
      }
    }
  }
  return intersec;
}

vec3 findColorBBOX (intersection hit, ray r, int count, bbox b) {
  // First look at point light source
  RGBQUAD toReturn;
  vec3 color;
  vec3 color1;
  vec3 color2;
  vec3 ambient;
  vec3 diffuse;
  vec3 specular;
  vec3 emission; // = vec3(emission[0], emission[1], emission[2]);
  float shiny;
  vec3 eyedirn = normalize (r.direction());
  if (pntExits == 1 && hit.getN() == 1) {
    vec3 pointcoords = vec3(point[0], point[1], point[2]);
    vec3 dir = normalize (hit.getPosition() - pointcoords);
    vec3 lightcolor = vec3(point[3], point[4], point[5]);
    vec3 normal;
    vec3 epsilon = dir * -0.001f;

    if (hit.getShape() == 1) { // Triangle
      diffuse = vec3(hit.getTriangle().diffuse[0], hit.getTriangle().diffuse[1], hit.getTriangle().diffuse[2]);
      specular = vec3(hit.getTriangle().specular[0], hit.getTriangle().specular[1], hit.getTriangle().specular[2]);
      shiny = hit.getTriangle().shininess;
      normal = hit.getNormal();
    } else if (hit.getShape() == 2) { // Sphere
      diffuse = vec3(hit.getSphere().diffuse[0], hit.getSphere().diffuse[1], hit.getSphere().diffuse[2]);
      specular = vec3(hit.getSphere().specular[0], hit.getSphere().specular[1], hit.getSphere().specular[2]);
      shiny = hit.getSphere().shininess;
      normal = hit.getNormal();
    }
    vec3 halfvec = normalize (dir + eyedirn);
    float distancetolight = sqrt(pow(hit.getPosition().x - pointcoords.x, 2) + pow(hit.getPosition().y - pointcoords.y, 2) + pow(hit.getPosition().z - pointcoords.z, 2));
    float attencalcdenom = attenuation[0] + (attenuation[1] * distancetolight) + (attenuation[2] * pow(distancetolight, 2));
    vec3 attenlight = lightcolor / attencalcdenom;
    if (attenExists != 1) {
      attenlight = lightcolor;
      //cout << "attennot working" << endl;
    }
    float nDotL = dot(normal, dir)  ;      
    // vec3 lambert = diffuse * attenlight * max (nDotL, 0.0f) ;  
    vec3 lambert = diffuse * attenlight * max (nDotL, 0.0f) ; 
    float nDotH = dot(normal, halfvec) ; 
    // vec3 phong = specular * attenlight * pow (max(nDotH, 0.0f), shiny) ; 
    vec3 phong = specular * attenlight * pow (max(nDotH, 0.0f), shiny) ; 
    color1 = lambert + phong;
    // ra ray from light source to point
    ray ra = ray(hit.getPosition() + epsilon, normalize(dir * -1.0f));
    // intersection obj = Intersect(ra);
    intersection obj = IntersectBBOX(ra, b);
    if (obj.getN() && obj.getDistance() < distancetolight) {
      color1 = vec3(0.0f, 0.0f, 0.0f);
    } 
  }
    //color += color1;
  if (dirExits == 1 && hit.getN() == 1) {
    vec3 pointcoords = vec3(directional[0], directional[1], directional[2]);
    vec3 dir = normalize (-1.0f * pointcoords);
    vec3 lightcolor = vec3(directional[3], directional[4], directional[5]);
    vec3 normal;
    vec3 epsilon = dir * -0.001f;

    if (hit.getShape() == 1) { // Triangle
      diffuse = vec3(hit.getTriangle().diffuse[0], hit.getTriangle().diffuse[1], hit.getTriangle().diffuse[2]);
      specular = vec3(hit.getTriangle().specular[0], hit.getTriangle().specular[1], hit.getTriangle().specular[2]);
      shiny = hit.getTriangle().shininess;
      normal = hit.getNormal();
    } else if (hit.getShape() == 2) { // Sphere
      diffuse = vec3(hit.getSphere().diffuse[0], hit.getSphere().diffuse[1], hit.getSphere().diffuse[2]);
      specular = vec3(hit.getSphere().specular[0], hit.getSphere().specular[1], hit.getSphere().specular[2]);
      shiny = hit.getSphere().shininess;
      normal = hit.getNormal();
    }
    //float distancetolight = sqrt(pow(hit.getPosition().x - pointcoords.x, 2) + pow(hit.getPosition().y - pointcoords.y, 2) + pow(hit.getPosition().z - pointcoords.z, 2));
    vec3 halfvec = normalize (dir + eyedirn);
    float nDotL = dot(normal, dir)  ;      
    vec3 lambert = diffuse * lightcolor * max (nDotL, 0.0f) ;  
    float nDotH = dot(normal, halfvec) ; 
    vec3 phong = specular * lightcolor * pow (max(nDotH, 0.0f), shiny) ; 
    color2 = lambert + phong;
    // ra ray from light source to point
    ray ra = ray(hit.getPosition() + epsilon, normalize(dir * -1.0f));
    // intersection obj = Intersect(ra);
    intersection obj = IntersectBBOX(ra, b);
    if (obj.getN()) {
      color2 = vec3(0.0f, 0.0f, 0.0f);
    } 
  }
  color = color1 + color2;
  //cout << "color r: " << color.x << " color g: " << color.y << " color b: " << color.z << endl;
  return color;
  // cout << lambert[0] << " " << diffuse[1] << " " << diffuse[2] << endl;
  // toReturn.rgbRed = color[0] * 255;
  // toReturn.rgbGreen = color[1] * 255;
  // toReturn.rgbBlue = color[2] * 255;
  // return toReturn;
  
}

vec3 recursionBBOX(const ray& r, int depth, bbox b) {
  if (depth <= 0) {
    return vec3(0.0f); 
  }
  // intersection hit = Intersect(r);
  intersection hit = IntersectBBOX(r, b);
  if (hit.getN()) {
    vec3 directColor = findColorBBOX(hit, r, 0, b);
    vec3 reflectionColor = vec3(0.0f);
    vec3 specular;
    if (hit.getShape() == 1 || hit.getShape() == 2) { 
      vec3 reflectDir = reflect(r.direction(), hit.getNormal());
      vec3 epsilon = normalize(reflectDir) * 0.001f;
      ray reflectedRay(hit.getPosition() + epsilon, reflectDir);
      if (hit.getShape() == 1) {
        specular = vec3(hit.getTriangle().specular[0], hit.getTriangle().specular[1], hit.getTriangle().specular[2]);
      } else {
        specular = vec3(hit.getSphere().specular[0], hit.getSphere().specular[1], hit.getSphere().specular[2]);
      }
        reflectionColor = recursionBBOX(reflectedRay, depth - 1, b);
    }
    vec3 finalColor = directColor + reflectionColor * specular;


    return finalColor;
  } else {
    return vec3(0.0f, 0.0f, 0.0f);
  }
}






















intersection Intersect (const ray r) {
  float mindist = INFINITY;
  RGBQUAD color;
  intersection intersec;

  // Now checking spheres
  for (int i = 0; i < numspheres; i++) {
    //cout << "Should not run" << endl;
    vec3 cent = spheres[i].v;
    mat4 inv = inverse(spheres[i].trans);
    vec4 centtranstemp = inv * vec4(cent, 1.0f);
    vec4 rayoritemp = inv * vec4(r.origin(), 1.0f);
    vec4 raydirtemp = inv * vec4(r.direction(), 0.0f);
    vec3 centtrans = vec3(centtranstemp.x, centtranstemp.y, centtranstemp.z);
    vec3 rayori = vec3(rayoritemp.x, rayoritemp.y, rayoritemp.z);
    vec3 raydir = vec3(raydirtemp.x, raydirtemp.y, raydirtemp.z);
    ray raytrans = ray (rayori, raydir);
    // float a = dot(r.direction() , r.direction());
    // float b = 2 * dot(r.direction(), r.origin() - cent);
    // float c = dot(r.origin() - cent, r.origin() - cent) - (spheres[i].radius * spheres[i].radius);
    float a = dot(raydir , raydir);
    float b = 2 * dot(raydir, rayori - cent);
    float c = dot(rayori - cent, rayori - cent) - (spheres[i].radius * spheres[i].radius);
    float discriminant = (b * b) - (4 * a * c);
    //cout << "raydirection: " << " " << r.direction().x << " " << r.direction().y << " " <<r.direction().z << '\n';
    //cout << "a = " << a << " b = " << b << " c = " << c << endl;
    //cout << discriminant << endl;
    if (discriminant >= 0) {
      //cout << "potential" << endl;
      float t1 = (-b + sqrt(discriminant)) / (2 * a);
      float t2 = (-b - sqrt(discriminant)) / (2 * a);
      // Case both pos take the smaller of the two;
      if (t1 > 0 && t2 > 0) {
        //cout << "potential" << endl;
        float t = min (t1, t2);
        if (t < mindist) {
          mindist = t;
          vec3 point = rayori + (raydir * mindist);
          vec4 pointtemp = vec4(point, 1.0f);
          vec4 pointtemp2 = spheres[i].trans * pointtemp;
          vec3 pointposition = vec3(pointtemp2.x, pointtemp2.y, pointtemp2.z);
          // vec3 ntransnorm = r.origin() + (r.direction() * mindist);
          vec4 normtemp = inverseTranspose(spheres[i].trans) * vec4((point), 1.0f);
          vec3 norm = normalize (vec3(normtemp.x, normtemp.y, normtemp.z)) * -1.0f;
          // vec3 normpoint = r.origin() + (r.direction() * mindist) - cent;
          // vec4 normpointtrans = inv * vec4(normpoint, 1.0f);
          // vec3 norm = normalize (vec3(normpointtrans.x, normpointtrans.y, normpointtrans.z)) * -1.0f;
          intersec = intersection(pointposition, norm, mindist, 2);
          intersec.setSphere(spheres[i]);
        }
      // Case that one is pos
      } else if (t1 > 0 || t2 > 0) {
        //cout << "potential" << endl;
        float t = max (t1, t2);
        if (t < mindist) {
          mindist = t;
          vec3 point = rayori + (raydir * mindist);
          vec4 pointtemp = vec4(point, 1.0f);
          vec4 pointtemp2 = spheres[i].trans * pointtemp;
          vec3 pointposition = vec3(pointtemp2.x, pointtemp2.y, pointtemp2.z);
          // vec3 ntransnorm = r.origin() + (r.direction() * mindist);
          vec4 normtemp = inverseTranspose(spheres[i].trans) * vec4((point), 1.0f);
          vec3 norm = normalize (vec3(normtemp.x, normtemp.y, normtemp.z)) * -1.0f;
          // vec3 normpoint = r.origin() + (r.direction() * mindist) - cent;
          // vec4 normpointtrans = inv * vec4(normpoint, 1.0f);
          // vec3 norm = normalize (vec3(normpointtrans.x, normpointtrans.y, normpointtrans.z)) * -1.0f;
          intersec = intersection(pointposition, norm, mindist, 2);
          intersec.setSphere(spheres[i]);
        }
      }
      //cout << "t1: " << t1 << " t2: " << t2 << endl;
    }
  }
  // Now checking triangles
  
  for (int i = 0; i < numtriangles; i++) {
    //cout << "Should run twice" << endl;
    vec4 tempa = triangles[i].trans * vec4(triangles[i].v1, 1.0f);
    vec4 tempb = triangles[i].trans * vec4(triangles[i].v2, 1.0f);
    vec4 tempc = triangles[i].trans * vec4(triangles[i].v3, 1.0f);
    // vec3 a = triangles[i].v1;
    // vec3 b = triangles[i].v2;
    // vec3 c = triangles[i].v3;
    vec3 a = vec3(tempa.x, tempa.y, tempa.z);
    vec3 b = vec3(tempb.x, tempb.y, tempb.z);
    vec3 c = vec3(tempc.x, tempc.y, tempc.z);
    vec3 norm = normalize(cross(c-a, b-a));
    float denomcheck = dot(r.direction(), norm);
    if (denomcheck != 0) {
      float t = ((dot(a, norm) - dot(r.origin(), norm)) / dot(r.direction(), norm));
      vec3 pointp = r.origin() + (r.direction() * t);
      vec3 pminusa = pointp - a;
      vec3 bminusa = b - a;
      vec3 cminusa = c - a;
      float gammanumerator = (pminusa.y * bminusa.x) - (pminusa.x * bminusa.y);
      float gammadenominator = (bminusa.y * cminusa.x) - (cminusa.y * bminusa.x);
      float gamma = -1.0f * (gammanumerator / gammadenominator);
      float betapt1 = (pminusa.x / bminusa.x);
      float betapt2 = (gamma * cminusa.x) / bminusa.x;
      float beta = betapt1 - betapt2;
      float alpha = (pointp.x - (beta * b.x) - (gamma * c.x)) / a.x;
      float total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < mindist) || (t > 0 && mindist == -1)) {
          mindist = t;
          color.rgbRed = triangles[i].ambient[0] * 255;
          color.rgbGreen = triangles[i].ambient[1] * 255;
          color.rgbBlue = triangles[i].ambient[2] * 255;
          intersec = intersection((r.origin() + (r.direction() * mindist)), norm, mindist, 1);
          intersec.setTriangle(triangles[i]);
          //cout << "tri found" << endl;
        }
      } 
      gammanumerator = (pminusa.z * bminusa.y) - (pminusa.y * bminusa.z);
      gammadenominator = (bminusa.z * cminusa.y) - (cminusa.z * bminusa.y);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.y / bminusa.y);
      betapt2 = (gamma * cminusa.y) / bminusa.y;
      beta = betapt1 - betapt2;
      alpha = (pointp.y - (beta * b.y) - (gamma * c.y)) / a.y;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < mindist) || (t > 0 && mindist == -1)) {
          mindist = t;

          intersec = intersection((r.origin() + (r.direction() * mindist)), norm, mindist, 1);
          intersec.setTriangle(triangles[i]);
          //cout << "tri found" << endl;
        }
      }
      gammanumerator = (pminusa.x * bminusa.z) - (pminusa.z * bminusa.x);
      gammadenominator = (bminusa.x * cminusa.z) - (cminusa.x * bminusa.z);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.z / bminusa.z);
      betapt2 = (gamma * cminusa.z) / bminusa.z;
      beta = betapt1 - betapt2;
      alpha = (pointp.z - (beta * b.z) - (gamma * c.z)) / a.z;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < mindist) || (t > 0 && mindist == -1)) {
          mindist = t;

          intersec = intersection((r.origin() + (r.direction() * mindist)), norm, mindist, 1);
          intersec.setTriangle(triangles[i]);
          //cout << "tri found" << endl;
        }
      }
      gammanumerator = (pminusa.x * bminusa.y) - (pminusa.y * bminusa.x);
      gammadenominator = (bminusa.x * cminusa.y) - (cminusa.x * bminusa.y);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.y / bminusa.y);
      betapt2 = (gamma * cminusa.y) / bminusa.y;
      beta = betapt1 - betapt2;
      alpha = (pointp.y - (beta * b.y) - (gamma * c.y)) / a.y;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < mindist) || (t > 0 && mindist == -1)) {
          mindist = t;

          intersec = intersection((r.origin() + (r.direction() * mindist)), norm, mindist, 1);
          intersec.setTriangle(triangles[i]);
          //cout << "tri found" << endl;
        }
      } 
      gammanumerator = (pminusa.y * bminusa.z) - (pminusa.z * bminusa.y);
      gammadenominator = (bminusa.y * cminusa.z) - (cminusa.y * bminusa.z);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.z / bminusa.z);
      betapt2 = (gamma * cminusa.z) / bminusa.z;
      beta = betapt1 - betapt2;
      alpha = (pointp.z - (beta * b.z) - (gamma * c.z)) / a.z;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < mindist) || (t > 0 && mindist == -1)) {
          mindist = t;

          intersec = intersection((r.origin() + (r.direction() * mindist)), norm, mindist, 1);
          intersec.setTriangle(triangles[i]);
          //cout << "tri found" << endl;
        }
      }
      gammanumerator = (pminusa.z * bminusa.x) - (pminusa.x * bminusa.z);
      gammadenominator = (bminusa.z * cminusa.x) - (cminusa.z * bminusa.x);
      gamma = -1.0f * (gammanumerator / gammadenominator);
      betapt1 = (pminusa.x / bminusa.x);
      betapt2 = (gamma * cminusa.x) / bminusa.x;
      beta = betapt1 - betapt2;
      alpha = (pointp.x - (beta * b.x) - (gamma * c.x)) / a.x;
      total = beta + gamma; //+ alpha;
      if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && total <= 1) {
        if ((t > 0 && t < mindist) || (t > 0 && mindist == -1)) {
          mindist = t;

          intersec = intersection((r.origin() + (r.direction() * mindist)), norm, mindist, 1);
          intersec.setTriangle(triangles[i]);
          //cout << "tri found" << endl;
        }
      }
    }
  }
  
  //cout << "mindist = " << mindist << endl;
  if (mindist == -1) {
    intersec = intersection();
  } 
  //cout << "Should return something" << endl;
  return intersec;
}

vec3 findColor (intersection hit, ray r, int count) {
  // First look at point light source
  RGBQUAD toReturn;
  vec3 color;
  vec3 color1;
  vec3 color2;
  vec3 ambient;
  vec3 diffuse;
  vec3 specular;
  vec3 emission; // = vec3(emission[0], emission[1], emission[2]);
  float shiny;
  vec3 eyedirn = normalize (r.direction());
  if (pntExits == 1 && hit.getN() == 1) {
    vec3 pointcoords = vec3(point[0], point[1], point[2]);
    vec3 dir = normalize (hit.getPosition() - pointcoords);
    vec3 lightcolor = vec3(point[3], point[4], point[5]);
    vec3 normal;
    vec3 epsilon = dir * -0.001f;

    if (hit.getShape() == 1) { // Triangle
      diffuse = vec3(hit.getTriangle().diffuse[0], hit.getTriangle().diffuse[1], hit.getTriangle().diffuse[2]);
      specular = vec3(hit.getTriangle().specular[0], hit.getTriangle().specular[1], hit.getTriangle().specular[2]);
      shiny = hit.getTriangle().shininess;
      normal = hit.getNormal();
    } else if (hit.getShape() == 2) { // Sphere
      diffuse = vec3(hit.getSphere().diffuse[0], hit.getSphere().diffuse[1], hit.getSphere().diffuse[2]);
      specular = vec3(hit.getSphere().specular[0], hit.getSphere().specular[1], hit.getSphere().specular[2]);
      shiny = hit.getSphere().shininess;
      normal = hit.getNormal();
    }
    vec3 halfvec = normalize (dir + eyedirn);
    float distancetolight = sqrt(pow(hit.getPosition().x - pointcoords.x, 2) + pow(hit.getPosition().y - pointcoords.y, 2) + pow(hit.getPosition().z - pointcoords.z, 2));
    float attencalcdenom = attenuation[0] + (attenuation[1] * distancetolight) + (attenuation[2] * pow(distancetolight, 2));
    vec3 attenlight = lightcolor / attencalcdenom;
    if (attenExists != 1) {
      attenlight = lightcolor;
      //cout << "attennot working" << endl;
    }
    float nDotL = dot(normal, dir)  ;      
    // vec3 lambert = diffuse * attenlight * max (nDotL, 0.0f) ;  
    vec3 lambert = diffuse * attenlight * max (nDotL, 0.0f) ; 
    float nDotH = dot(normal, halfvec) ; 
    // vec3 phong = specular * attenlight * pow (max(nDotH, 0.0f), shiny) ; 
    vec3 phong = specular * attenlight * pow (max(nDotH, 0.0f), shiny) ; 
    color1 = lambert + phong;
    // ra ray from light source to point
    ray ra = ray(hit.getPosition() + epsilon, normalize(dir * -1.0f));
    // intersection obj = Intersect(ra);
    intersection obj = Intersect(ra);
    if (obj.getN() && obj.getDistance() < distancetolight) {
      color1 = vec3(0.0f, 0.0f, 0.0f);
    } 
  }
    //color += color1;
  if (dirExits == 1 && hit.getN() == 1) {
    vec3 pointcoords = vec3(directional[0], directional[1], directional[2]);
    vec3 dir = normalize (-1.0f * pointcoords);
    vec3 lightcolor = vec3(directional[3], directional[4], directional[5]);
    vec3 normal;
    vec3 epsilon = dir * -0.001f;

    if (hit.getShape() == 1) { // Triangle
      diffuse = vec3(hit.getTriangle().diffuse[0], hit.getTriangle().diffuse[1], hit.getTriangle().diffuse[2]);
      specular = vec3(hit.getTriangle().specular[0], hit.getTriangle().specular[1], hit.getTriangle().specular[2]);
      shiny = hit.getTriangle().shininess;
      normal = hit.getNormal();
    } else if (hit.getShape() == 2) { // Sphere
      diffuse = vec3(hit.getSphere().diffuse[0], hit.getSphere().diffuse[1], hit.getSphere().diffuse[2]);
      specular = vec3(hit.getSphere().specular[0], hit.getSphere().specular[1], hit.getSphere().specular[2]);
      shiny = hit.getSphere().shininess;
      normal = hit.getNormal();
    }
    //float distancetolight = sqrt(pow(hit.getPosition().x - pointcoords.x, 2) + pow(hit.getPosition().y - pointcoords.y, 2) + pow(hit.getPosition().z - pointcoords.z, 2));
    vec3 halfvec = normalize (dir + eyedirn);
    float nDotL = dot(normal, dir)  ;      
    vec3 lambert = diffuse * lightcolor * max (nDotL, 0.0f) ;  
    float nDotH = dot(normal, halfvec) ; 
    vec3 phong = specular * lightcolor * pow (max(nDotH, 0.0f), shiny) ; 
    color2 = lambert + phong;
    // ra ray from light source to point
    ray ra = ray(hit.getPosition() + epsilon, normalize(dir * -1.0f));
    // intersection obj = Intersect(ra);
    intersection obj = Intersect(ra);
    if (obj.getN()) {
      color2 = vec3(0.0f, 0.0f, 0.0f);
    } 
  }
  color = color1 + color2;
  //cout << "color r: " << color.x << " color g: " << color.y << " color b: " << color.z << endl;
  return color;
  // cout << lambert[0] << " " << diffuse[1] << " " << diffuse[2] << endl;
  // toReturn.rgbRed = color[0] * 255;
  // toReturn.rgbGreen = color[1] * 255;
  // toReturn.rgbBlue = color[2] * 255;
  // return toReturn;
  
}

vec3 recursion(const ray& r, int depth) {
  if (depth <= 0) {
    return vec3(0.0f); 
  }
  // intersection hit = Intersect(r);
  intersection hit = Intersect(r);
  if (hit.getN()) {
    vec3 directColor = findColor(hit, r, 0);
    vec3 reflectionColor = vec3(0.0f);
    vec3 specular;
    if (hit.getShape() == 1 || hit.getShape() == 2) { 
      vec3 reflectDir = reflect(r.direction(), hit.getNormal());
      vec3 epsilon = normalize(reflectDir) * 0.001f;
      ray reflectedRay(hit.getPosition() + epsilon, reflectDir);
      if (hit.getShape() == 1) {
        specular = vec3(hit.getTriangle().specular[0], hit.getTriangle().specular[1], hit.getTriangle().specular[2]);
      } else {
        specular = vec3(hit.getSphere().specular[0], hit.getSphere().specular[1], hit.getSphere().specular[2]);
      }
        reflectionColor = recursion(reflectedRay, depth - 1);
    }
    vec3 finalColor = directColor + reflectionColor * specular;


    return finalColor;
  } else {
    return vec3(0.0f, 0.0f, 0.0f);
  }
}

ray RayThruPixel (const vec3& thisw, const vec3& u, const vec3& v, int i, int j) {
  
  float alpha = tan((xfov) / 2.0f) * (((float)j + 0.5 - ((float)width / 2.0f)) / ((float)width / 2.0f)); //horizontal pixels
  float beta = tan((yfov) / 2.0f) * ((((float)height / 2.0f) - (float)i + 0.5) / ((float)height / 2.0f)); //vertical pixels
  

  ray r = ray(eye, normalize((alpha * u) + (beta * v) - thisw));
  //cout << tan((xfov)/2.0f) << '\n';
  //cout << ((((float)height / 2.0f) - (float)i) / ((float)height / 2.0f)) << '\n';
  //cout << j<<endl;
  //cout << "alpha = " << alpha << " beta = " << beta << endl;
  return r;
  //return vec3(0,0,0);//return ray;
}



// void Raytrace (const vec3& eye, const vec3& center, const vec3& up, RGBQUAD im[480][640]) {
//   //Image image = new image (width, height);
//   for (int i = 0; i < height; i++) {
//     for (int j = 0; j < width; j++) {
//       ray r = RayThruPixel (eye, center, up, i, j);
//       im[i][j] = Intersect (r, eye);
//       //cout << "ran" << endl;
//       //image[i][j] = FindColor (hit);
//     }
//   }
//   cout << "finished" << endl;
//   //return image;
// }


int main(int argc, char* argv[]) {
  readfile(argv[1]);
  //cout << numtriangles << endl;
  FreeImage_Initialise();
  //RGBQUAD im[500][500];
  //Raytrace(eye, center, up, im);

  FIBITMAP *bitmap = FreeImage_Allocate(width, height, 24);

  //cout << "eye: " << eye.x << " " << eye.y << " " << eye.z << '\n';
  //cout << "center: " << center.x << " " << center.y << " " << center.z << '\n';
  //cout << "up: " << up.x << " " << up.y << " " << up.z << '\n';

  vec3 thisw = normalize(eye - center); // named thisw to not mix up with w for width
  vec3 temp = cross(up, thisw);
  vec3 u = normalize(temp);
  vec3 v = cross(thisw, u);
    //RGBQUAD color ;
  // for (int j =0; j <width; j++) {
  object obs[pilesize];
  int obssize = 0;
  for (int i = 0; i < numtriangles; i++) {
    object temp = object(1);
    temp.setTriangle(triangles[i]);
    obs[i] = temp;
    obssize++;
  }
  for (int i = numtriangles; i < numtriangles + numspheres; i++) {
    object temp = object(2);
    temp.setSphere(spheres[i - numtriangles]);
    obs[i] = temp;
    obssize++;
  }
  bbox b = bboxcreate(obs, obssize, 0);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      ray r = RayThruPixel(thisw, u, v, i, j);
      // float p1 = r.direction().x;
      // p1 = (p1 * 0.5f + 0.5f) * 255.0F;
      // RGBQUAD color;
      // color.rgbBlue = p1;
      // FreeImage_SetPixelColor(bitmap,j,i,&color);
      // intersection hit = Intersect(r);
      intersection hit = IntersectBBOX(r, b);
      RGBQUAD a;
      if (hit.getN()) {
        vec3 color = recursionBBOX(r, 4, b);
        vec3 ambient;
        vec3 emission;
        if (hit.getShape() == 1) { // Triangle
          ambient = vec3(hit.getTriangle().ambient[0], hit.getTriangle().ambient[1], hit.getTriangle().ambient[2]);
          emission = vec3(hit.getTriangle().emission[0], hit.getTriangle().emission[1], hit.getTriangle().emission[2]);
        } else if (hit.getShape() == 2) { // Sphere
          ambient = vec3(hit.getSphere().ambient[0], hit.getSphere().ambient[1], hit.getSphere().ambient[2]);
          emission = vec3(hit.getSphere().emission[0], hit.getSphere().emission[1], hit.getSphere().emission[2]);
        }
        color = color + ambient + emission;
        if (color.x > 1.0f || color.y > 1.0f || color.z > 1.0f) {
          //cout << color.x << " " << color.y << " " << color.z << endl;
        }
        a.rgbRed = color.x * 255.0f;
        a.rgbGreen = color.y * 255.0f;
        a.rgbBlue = color.z * 255.0f;
      } else {
        a.rgbRed = 0;
        a.rgbGreen = 0;
        a.rgbBlue = 0;
      }
      //RGBQUAD a = findColor(hit, r, 0);
      
      //FreeImage_SetPixelColor(bitmap, height - i, j, &a);
      FreeImage_SetPixelColor(bitmap, j, height - i, &a);
      //for (int k = 0; k < numtriangles; k++) {
      //FreeImage_SetPixelColor(bitmap, height - i, j, &a);
      //}
      //FreeImage_SetPixelColor(bitmap, i, j, &im[j][i]);
    }
  }


    //FIBITMAP img = FreeImage_ConvertFromRawBits(pixels, w, h, w 3, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);


    //FreeImage_Save(FIF_PNG, img, "img.png", 0);
  // int pix = w * h;
  // BYTE *pixels = new BYTE[3*pix];	

  FreeImage_Save(FIF_PNG, bitmap, argv[2], 0);
  

  FreeImage_DeInitialise();
  return 0;
}
