// Readfile definitions 
#include <sstream>
#include <stack>
using namespace std; 
//#include <C:\Users\mmizu\OneDrive\Desktop\School\CSE167\hw4-windows\hw2-windows\hw2-windows\packages\glm.0.9.7.1\build\native\include\glm\glm.hpp>
//#include <C:\Users\mmizu\OneDrive\Desktop\School\CSE167\hw4-windows\hw2-windows\hw2-windows\packages\glm.0.9.7.1\build\native\include\glm\gtc\matrix_transform.hpp>

void matransform (stack<mat4> &transfstack, float * values) ;
void rightmultiply (const mat4 & M, stack<mat4> &transfstack) ;
bool readvals (stringstream &s, const int numvals, float * values) ;
void readfile (const char * filename) ;
