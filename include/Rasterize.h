#ifndef _RASTERIZE_H_
#define _RASTERIZE_H_


#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

#include "Model_Transformations.h"
#include "View_Transformations.h"


using namespace tinyxml2;
using namespace std;


bool isVisible(double den, double num, double &tE, double &tL);

void interpolatePoint(Vec4 &point, const Vec4 &start, double t, double dx, double dy, double dz);

void interpolateColor(Color &color, const Color &start, const Color &diff_color, double t);

vector<Vec4> processTrianglePoints(const Triangle &triangle, const Matrix4 &matrix4, const std::vector<Vec3 *> &vertices);

vector<Color> getTriangleColors(const Triangle &triangle, const std::vector<Vec3 *> &vertices, const Scene *scene);

void rasterizeLine(Vec4& first_point, Vec4& second_point, Color& first_color, Color& second_color, Matrix4 viewport_matrix4, std::vector<std::vector<Color> >& image);

void rasterizeSolid(vector<Vec4> pVector, vector<Color> color_Vector, Matrix4 pMatrix4, Camera *camera, std::vector<std::vector<Color>> &image , Scene*scene);


#endif