#ifndef _MODEL_TRANSFORMATION_H_
#define _MODEL_TRANSFORMATION_H_

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

#include "View_Transformations.h"

using namespace tinyxml2;
using namespace std;


#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"





Matrix4 getModelMat(Mesh *mesh, std::vector<Scaling *>  scalings , std::vector<Rotation *> rotations , std::vector<Translation *> translations );






#endif