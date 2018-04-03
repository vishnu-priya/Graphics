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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <iostream>

using namespace std;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}
MGLpixel colorGiven[3];

class Matrix{
public:
  MGLfloat tm[4][4];
  Matrix():tm{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
  {}
  Matrix(MGLfloat X, MGLfloat Y, MGLfloat Z, MGLfloat w){
    tm[0][0] = X; tm[1][1] = Y; tm[2][2] = Z; tm[3][3] = w;
  }
  Matrix& operator = (const Matrix& m){
		if (this != &m){
			for (int i=0; i<4; i++)
				for (int j=0; j<4; j++)
					tm[i][j] = m.tm[i][j];
		}
    return *this;
	}
  Matrix& mult(const Matrix& m, const Matrix& x){
    if(this != &m){
      for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
        for(int k=0; k<4; k++){
          tm[i][j] += m.tm[i][k] * x.tm[k][j];
      }
    }
    return *this;
  }
};

class vertex {
public:
  MGLfloat x, y, z, w;
  MGLpixel color[3];
  vertex (MGLfloat a)
  { x = a; y=0; z=0;}
  vertex (MGLfloat a, MGLfloat b, MGLpixel* c)
  { x= a; y= b; z= 0; color[0]= c[0]; color[1]= c[1]; color[2]= c[2];}
  vertex (MGLfloat a, MGLfloat b, MGLfloat d, MGLpixel* c)
  { x= a; y= b; z= d; color[0]= c[0]; color[1]= c[1]; color[2]= c[2];}
  MGLfloat area(vertex a, vertex b, vertex c){
    return ((b.y-c.y)*a.x + (c.x-b.x)*a.y + b.x*c.y - c.x*b.y)/2;
  }
};

class pixel{
public:
  MGLfloat x, y;
  MGLpixel color;
  MGLfloat zbuffer;
  pixel (MGLfloat X, MGLfloat Y, MGLpixel c, MGLfloat z)
  {x=X; y=Y; color=c; zbuffer=z;}
};
struct shapes{
  MGLsize N;
};
std::vector<shapes> nshapes;
int mgl_MatrixMode;
Matrix ortho;
Matrix currentMatrix = Matrix(1, 1, 1, 1);
std::vector<Matrix> currentMatrixStackModel, currentMatrixStackProjection;
std::vector<pixel> frameBuffer;
std::vector<vertex> v;

void drawTriangle(vertex& a,vertex& b,vertex& c)
{
  int x_max = floor(std::max({a.x, b.x, c.x}));
  int x_min = ceil(std::min({a.x, b.x, c.x}));
  int y_max = floor(std::max({a.y, b.y, c.y}));
  int y_min = ceil(std::min({a.y, b.y, c.y}));
  MGLfloat areaT = a.area(a, b, c);
  for(int i=x_min; i<=x_max; i++)
  for(int j=y_min; j<=y_max; j++){
    vertex P = vertex(MGLfloat(i+0.5), MGLfloat(j+0.5), colorGiven);
    MGLfloat alpha = a.area(P, b, c)/areaT;
    MGLfloat beta = b.area(P, c, a)/areaT;
    MGLfloat gamma = c.area(P, a, b)/areaT;
    if( alpha >=0 && alpha<=1 &&beta >=0 && beta <=1 && gamma >=0 && gamma <= 1){
      MGLfloat z = alpha*a.z+beta*b.z+gamma*c.z;
      MGLpixel color[3];
      for(int d=0; d<3; d++){
        color[d] = a.color[d]*alpha+b.color[d]*beta+c.color[d]*gamma;
      }
      MGLpixel col = Make_Pixel(color[0], color[1], color[2]);
      frameBuffer.push_back(pixel(i, j, col, z));
    }
  }
}
bool sortByZ(pixel a, pixel b) { return a.zbuffer > b.zbuffer; }

/**
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
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
  for(int i =0; i< height*width; i++)
    data[i] = Make_Pixel(1, 1, 1);
  for(int s=0; s<nshapes.size(); s++){
    int preset = s*nshapes[s-1].N;
    for(int i=preset; i< preset+nshapes[s].N; i++){
      v[i].x = width*(v[i].x)/2 +width/2-0.5;
    }
    for(int i=preset; i< preset+nshapes[s].N; i++){
      v[i].y = height*(v[i].y)/2 +height/2-0.5;
    }
    if(nshapes[s].N==3)
      drawTriangle(v[preset+0], v[preset+1], v[preset+2]);
    else if(nshapes[s].N==4){
      drawTriangle(v[preset+0], v[preset+1], v[preset+2]);
      drawTriangle(v[preset+0], v[preset+3], v[preset+2]);
    }
  }
  sort(frameBuffer.begin(), frameBuffer.end(), sortByZ);
	 for (int i = 0; i < frameBuffer.size(); ++i){
		int x = frameBuffer[i].x;
		int y = frameBuffer[i].y;
		MGLpixel color = frameBuffer[i].color;
		data[y*width + x] = color;
	}
  ortho = Matrix(0, 0, 0, 0);
  currentMatrix = Matrix(0, 0, 0, 0);
  currentMatrixStackModel.clear(); currentMatrixStackProjection.clear(); frameBuffer.clear();
  nshapes.clear(); v.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
  int sum = 0;
  for(int i=0; i<4; i++)
  for(int j=0; j<4; j++)
    sum += ortho.tm[i][j];
  if(sum==0){
    mglOrtho(-1, 1, -1, 1, -1, 1);
  }
  shapes s;
  switch(mode){
    case MGL_TRIANGLES:
      s.N = 3;
      nshapes.push_back(s);
      break;
    case MGL_QUADS:
      s.N=4;
      nshapes.push_back(s);
      break;
  }
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  Matrix transform;
  transform = transform.mult(ortho, currentMatrix);
  int count = 0;
  for(int i=0; i<nshapes.size(); i++){
    count += nshapes[i].N;
  }
  int preset = 0;
  if(nshapes.size()>1){
    preset = nshapes[nshapes.size()-1].N+1;
  }
  if(count != v.size()){
    for(int i=count; i<v.size(); i+=nshapes.back().N){
      shapes s;
      s.N = nshapes.back().N;
      nshapes.push_back(s);
    }
  }
  for(int i=preset; i<v.size(); i++){
    v[i].w = v[i].x*transform.tm[3][0]+v[i].y*transform.tm[3][1]+v[i].z*transform.tm[3][2]+transform.tm[3][3];
  }
  for(int i=preset; i< v.size(); i++){
    v[i].x = (v[i].x*transform.tm[0][0]+v[i].y*transform.tm[0][1]+v[i].z*transform.tm[0][2]+transform.tm[0][3])/v[i].w;
  }
  for(int i=preset; i< v.size(); i++){
    v[i].y = (v[i].x*transform.tm[1][0]+v[i].y*transform.tm[1][1]+ v[i].z*transform.tm[1][2] + transform.tm[1][3])/v[i].w;
  }
  for(int i=preset; i< v.size(); i++){
    v[i].z = (v[i].x*transform.tm[2][0]+ v[i].y*transform.tm[2][1]+v[i].z*transform.tm[2][2] + transform.tm[2][3])/v[i].w;
  }
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
  v.push_back(vertex(x, y, colorGiven));
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  v.push_back(vertex(x, y, z, colorGiven));
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  if (mode == MGL_MODELVIEW || mode == MGL_PROJECTION)
		mgl_MatrixMode = mode;
	else
		mgl_MatrixMode = -1;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
  if (mgl_MatrixMode == MGL_MODELVIEW && !currentMatrixStackModel.empty())
		currentMatrixStackModel.push_back(currentMatrixStackModel.back());
	else if (mgl_MatrixMode == MGL_PROJECTION && !currentMatrixStackProjection.empty())
		currentMatrixStackProjection.push_back(currentMatrixStackProjection.back());
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
  if (mgl_MatrixMode == MGL_MODELVIEW && !currentMatrixStackModel.empty()){
		currentMatrixStackModel.pop_back();
    currentMatrix = currentMatrixStackModel.back();
	}else if (mgl_MatrixMode == MGL_PROJECTION && !currentMatrixStackProjection.empty()){
		currentMatrixStackProjection.pop_back();
    currentMatrix = currentMatrixStackProjection.back();
  }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  Matrix temp;
  for(int i=0; i<4; i++)
  for(int j=0; j<4; j++){
    if(i==j){
      temp.tm[i][j] = 1;
    }
  }
  if (mgl_MatrixMode == MGL_MODELVIEW){
    if(!currentMatrixStackModel.empty())
      currentMatrixStackModel.pop_back();
    currentMatrixStackModel.push_back(temp);
	}else if (mgl_MatrixMode == MGL_PROJECTION){
    if(!currentMatrixStackProjection.empty())
      currentMatrixStackProjection.pop_back();
    currentMatrixStackProjection.push_back(temp);
  }
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
  Matrix mult;
  int count = 0;
  for(int i=0; i<4;i++)
  for(int j=0; j<4; j++)
    mult.tm[j][i] = matrix[count++];
  currentMatrix = mult;
  if (mgl_MatrixMode == MGL_MODELVIEW)
		currentMatrixStackModel.push_back(mult);
	else if (mgl_MatrixMode == MGL_PROJECTION)
		currentMatrixStackProjection.push_back(mult);
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
  Matrix m, mult;
  if (mgl_MatrixMode == MGL_MODELVIEW && !currentMatrixStackModel.empty()){
		m = currentMatrixStackModel.back();
    currentMatrixStackModel.pop_back();
	}else if (mgl_MatrixMode == MGL_PROJECTION && !currentMatrixStackProjection.empty()){
		m = currentMatrixStackProjection.back();
    currentMatrixStackProjection.pop_back();
  }
  int count = 0;
  for(int i=0; i<4; i++){
    count = 0;
    for(int j=0; j<4; j++)
    for(int k=0; k<4; k++){
      mult.tm[i][j] += m.tm[i][k] * matrix[count++];
    }
  }
  currentMatrix = mult;
  if (mgl_MatrixMode == MGL_MODELVIEW)
		currentMatrixStackModel.push_back(mult);
	else if (mgl_MatrixMode == MGL_PROJECTION)
		currentMatrixStackProjection.push_back(mult);
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
  MGLfloat m[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1};
  mglMultMatrix(m);
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
  MGLfloat c = cos(angle* M_PI / 180);
  MGLfloat s = sin(angle* M_PI / 180);
  MGLfloat norm = sqrt(x*x+y*y+z*z);
  x /= norm; y /= norm; z /= norm;
  Matrix rotate;
  rotate.tm[0][0] = x*x*(1-c)+c;
  rotate.tm[0][1] = x*y*(1-c)-z*s;
  rotate.tm[0][2] = x*z*(1-c)+y*s;
  rotate.tm[1][0] = x*y*(1-c)+z*s;
  rotate.tm[1][1] = y*y*(1-c)+c;
  rotate.tm[1][2] = y*z*(1-c)-x*s;
  rotate.tm[2][0] = z*x*(1-c)-y*s;
  rotate.tm[2][1] = z*y*(1-c)+x*s;
  rotate.tm[2][2] = z*z*(1-c)+c;
  MGLfloat m[16] = {rotate.tm[0][0], rotate.tm[1][0], rotate.tm[2][0], 0, rotate.tm[0][1], rotate.tm[1][1], rotate.tm[2][1], 0, rotate.tm[0][2], rotate.tm[1][2], rotate.tm[2][2], 0, 0, 0, 0, 1};
  mglMultMatrix(m);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  MGLfloat m[16] = {x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1};
  mglMultMatrix(m);
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
  Matrix frustum;
  frustum.tm[0][0] = 2*near/(right-left);
  frustum.tm[1][1] = 2*near/(top-bottom);
  frustum.tm[0][2] = (right+left)/(right-left);
  frustum.tm[1][2] = (top+bottom)/(top-bottom);
  frustum.tm[2][2] = (far+near)/(near-far);
  frustum.tm[2][3] = (2*far*near)/(near-far);
  frustum.tm[3][2] = -1;
  MGLfloat m[16] = {frustum.tm[0][0], 0, 0, 0, 0, frustum.tm[1][1], 0, 0, frustum.tm[0][2], frustum.tm[1][2], frustum.tm[2][2], -1, 0, 0, frustum.tm[3][2], 0};
  mglMultMatrix(m);
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
  ortho.tm[0][0] = 2/(right-left);
  ortho.tm[1][1] = 2/(top-bottom);
  ortho.tm[2][2] = 2/(near-far);
  ortho.tm[3][3] = 1;
  ortho.tm[0][3] = (right+left)/(left-right);
  ortho.tm[1][3] = (top+bottom)/(bottom-top);
  ortho.tm[2][3] = (far+near)/(near-far);
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  colorGiven[0]= red*255; colorGiven[1]=green*255; colorGiven[2]=blue*255;
}
