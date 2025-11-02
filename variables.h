#ifdef MAINPROGRAM 
#define EXTERN 
#else
#define EXTERN extern 
#endif

EXTERN int amount; // The amount of rotation for each arrow press
EXTERN vec3 eye; // The (regularly updated) vector coordinates of the eye 
EXTERN vec3 up;  // The (regularly updated) vector coordinates of the up 

#ifdef MAINPROGRAM 
vec3 eyeinit(0.0,0.0,5.0) ; // Initial eye position, also for resets
vec3 upinit(0.0,1.0,0.0) ; // Initial up position, also for resets
vec3 center(0.0,0.0,0.0) ; // Center look at point 
int amountinit = 5;
float w = 500, h = 500 ; // width and height 
float fovy = 90.0 ; // For field of view
//float fovx = 90.0 ; 
#else 
EXTERN vec3 eyeinit ; 
EXTERN vec3 upinit ; 
EXTERN vec3 center ; 
EXTERN int amountinit;
EXTERN float w, h ; 
EXTERN float fovy ; 
//EXTERN float fovx ; 
#endif 

static enum {view, translate, scale} transop ; // which operation to transform 
//enum shape {triangle, sphere} ;
EXTERN float sx, sy ; // the scale in x and y 
EXTERN float tx, ty ; // the translation in x and y
EXTERN float width;
EXTERN float height;
EXTERN float yfov;
EXTERN float xfov;
EXTERN int depth;
//EXTERN float foxx;

// Lighting parameter array, similar to that in the fragment shader

// Materials (read from file) 
// With multiple objects, these are colors for each.
EXTERN int dirExits;
EXTERN float directional[6] ; 
EXTERN int pntExits;
EXTERN float point[6] ; 
EXTERN float ambient[3] ;

EXTERN float diffuse[3] ; 
EXTERN float specular[3] ; 
EXTERN float shininess ;
EXTERN float emission[3] ;
EXTERN float attenuation[3] ;
EXTERN int attenExists ;

const int pilesize = 150100;
EXTERN int maxverts;
EXTERN int numverts;
EXTERN struct vertex {
    vec3 v;
} vertexpile[pilesize]; 

EXTERN int maxvertnorms;
EXTERN int numvertnorms;
EXTERN struct vertexnormal {
    float x;
    float y;
    float z;
    float nx;
    float ny;
    float nz;
} vertexnormalpile[pilesize]; 

EXTERN int numtriangles;
EXTERN struct triangle {
  float ambient[3];
  float diffuse[3]; 
  float specular[3]; 
  float emission[3];
  float shininess;
  vec3 v1;
  vec3 v2;
  vec3 v3;
  mat4 trans;
} triangles[pilesize];

EXTERN int numtrianglenormals;
EXTERN struct trianglenormal {
  float ambient[3];
  float diffuse[3]; 
  float specular[3]; 
  float emission[3];
  float shininess;
  vertexnormal v1;
  vertexnormal v2;
  vertexnormal v3;
} trianglenormals[pilesize];

EXTERN int numspheres;
EXTERN struct sph {
  float ambient[3];
  float diffuse[3]; 
  float specular[3]; 
  float emission[3];
  float shininess;
  float radius;
  vec3 v;
  mat4 trans;
} spheres[pilesize];

// For multiple objects, read from a file.  
// const int maxobjects = 10 ; 
// EXTERN int numobjects ; 
// EXTERN struct object {
//   //shape type ; 
//   float size ;
//   float ambient[4] ; 
//   float diffuse[4] ; 
//   float specular[4] ;
//   float emission[4] ; 
//   float shininess ;
//   mat4 transform ; 
// } objects[maxobjects] ;
