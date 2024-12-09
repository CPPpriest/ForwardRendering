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


Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            vector<Color> rowOfColors;
            vector<double> rowOfDepths;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfColors.push_back(this->backgroundColor);
                rowOfDepths.push_back(std::numeric_limits<double>::infinity()); // Initialize depth buffer to maximum depth
            }

            this->image.push_back(rowOfColors);
            depthBuffer.push_back(rowOfDepths); // Initialize depth buffer
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
                depthBuffer[i][j] = std::numeric_limits<double>::infinity(); // Reset depth buffer
            }
        }
    }
}

int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}





void line_rasterization(Vec4 first_point, Vec4 second_point, Color first_color, Color second_color, Matrix4 viewport_matrix4, Camera *pCamera, std::vector<std::vector<Color> >& image) {
    Vec4 first_point_viewport = multiplyMatrixWithVec4(viewport_matrix4, first_point);
    Vec4 second_point_viewport = multiplyMatrixWithVec4(viewport_matrix4, second_point);

    // calculate the slope of the line
    double slope = (second_point_viewport.y - first_point_viewport.y) / (second_point_viewport.x - first_point_viewport.x);
    int x;
    int y;
    double x1_x0 = abs(second_point_viewport.x - first_point_viewport.x);
    double y1_y0 = abs(second_point_viewport.y - first_point_viewport.y);
    if(slope >= 0){
        x = min(first_point_viewport.x, second_point_viewport.x);
        y = min(first_point_viewport.y, second_point_viewport.y);
    }
    else if(slope < -1){
        x = max(first_point_viewport.x, second_point_viewport.x);
        y = min(first_point_viewport.y, second_point_viewport.y);
    }
    else {
        x = min(first_point_viewport.x, second_point_viewport.x);
        y = max(first_point_viewport.y, second_point_viewport.y);
    }
    if(slope < 1 && slope >= 0){
        double d = 2 * y1_y0 - x1_x0;
        double dE = 2 * y1_y0;
        double dNE = 2 * (y1_y0 - x1_x0);
        while(x <= max(first_point_viewport.x, second_point_viewport.x)){
 // ------------------------------- //
            image[x][y].r = round((first_color.r * abs(x - second_point_viewport.x) + second_color.r * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].g = round((first_color.g * abs(x - second_point_viewport.x) + second_color.g * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].b = round((first_color.b * abs(x - second_point_viewport.x) + second_color.b * abs(x - first_point_viewport.x)) / x1_x0);
 // ------------------------------- //            
            if(d <= 0){
                d += dE;
                x++;
            }
            else{
                d += dNE;
                x++;
                y++;
            }
        }
    }
    else if(slope >= 1){
        double d = 2 * x1_x0 - y1_y0;
        double dE = 2 * x1_x0;
        double dNE = 2 * (x1_x0 - y1_y0);
        while(y <= max(first_point_viewport.y, second_point_viewport.y)){
 // ------------------------------- //
            image[x][y].r = round((first_color.r * abs(y - second_point_viewport.y) + second_color.r * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].g = round((first_color.g * abs(y - second_point_viewport.y) + second_color.g * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].b = round((first_color.b * abs(y - second_point_viewport.y) + second_color.b * abs(y - first_point_viewport.y)) / y1_y0);
 // ------------------------------- //
            if(d <= 0){
                d += dE;
                y++;
            }
            else{
                d += dNE;
                x++;
                y++;
            }
        }
    }
    else if(slope < -1){
        double d = 2 * x1_x0 - y1_y0;
        double dE = 2 * x1_x0;
        double dNE = 2 * (x1_x0 - y1_y0);
        while(y < max(first_point_viewport.y, second_point_viewport.y)){
 // ------------------------------- //
            image[x][y].r = round((first_color.r * abs(y - second_point_viewport.y) + second_color.r * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].g = round((first_color.g * abs(y - second_point_viewport.y) + second_color.g * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].b = round((first_color.b * abs(y - second_point_viewport.y) + second_color.b * abs(y - first_point_viewport.y)) / y1_y0);
 // ------------------------------- //            
            
            if(d <= 0){
                d += dE;
                y++;
            }
            else{
                d += dNE;
                x--;
                y++;
            }
        }
    }
    else{
        double d = 2 * y1_y0 - x1_x0;
        double dE = 2 * y1_y0;
        double dNE = 2 * (y1_y0 - x1_x0);
        while(x < max(first_point_viewport.x, second_point_viewport.x)){
            image[x][y].r = round((first_color.r * abs(x - second_point_viewport.x) + second_color.r * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].g = round((first_color.g * abs(x - second_point_viewport.x) + second_color.g * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].b = round((first_color.b * abs(x - second_point_viewport.x) + second_color.b * abs(x - first_point_viewport.x)) / x1_x0);
            if(d <= 0){
                d += dE;
                x++;
            }
            else{
                d += dNE;
                x++;
                y--;
            }
        }
    }
}

void rasterization(vector<Vec4> pVector, vector<Color> color_Vector, Matrix4 pMatrix4, Camera *camera, std::vector<std::vector<Color>> &image , Scene*scene)
{
    for (int i = 0; i < pVector.size(); ++i) {
        pVector[i] = multiplyVec4WithScalar(pVector[i], 1 / pVector[i].t);
        pVector[i] = multiplyMatrixWithVec4(pMatrix4, pVector[i]);
    }

    // Take the minimum and maximum x and y values
    int minX = pVector[0].x;
    int maxX = pVector[0].x;
    int minY = pVector[0].y;
    int maxY = pVector[0].y;
    for (int i = 1; i < pVector.size(); ++i) {
        if (pVector[i].x < minX) {
            minX = pVector[i].x;
        }
        if (pVector[i].x > maxX) {
            maxX = pVector[i].x;
        }
        if (pVector[i].y < minY) {
            minY = pVector[i].y;
        }
        if (pVector[i].y > maxY) {
            maxY = pVector[i].y;
        }
    }

    if (minX < 0) minX = 0;
    if (minX > camera->horRes - 1) minX = camera->horRes - 1;
    if (minY < 0) minY = 0;
    if (minY > camera->verRes - 1) minY = camera->verRes - 1;
    if (maxX < 0) maxX = 0;
    if (maxX > camera->horRes - 1) maxX = camera->horRes - 1;
    if (maxY < 0) maxY = 0;
    if (maxY > camera->verRes - 1) maxY = camera->verRes - 1;

    // Create Vec3 vector using the pVector by discarding t value
    vector<Vec3> vec3_Vector;
    for (int i = 0; i < pVector.size(); ++i) {
        Vec3 temp = Vec3(pVector[i].x, pVector[i].y, pVector[i].z, pVector[i].colorId);
        vec3_Vector.push_back(temp);
    }

    // Rasterization algorithm using barycentric coordinates
    for (int i = minX; i <= maxX; i++) {
        for (int j = minY; j <= maxY; j++) {
            // Get the barycentric coordinates of vec3_Vector[0], vec3_Vector[1], vec3_Vector[2]
            
            
// ------------------------------- //
            double alpha = ((vec3_Vector[1].y - vec3_Vector[2].y) * i + (vec3_Vector[2].x - vec3_Vector[1].x) * j + vec3_Vector[1].x * vec3_Vector[2].y - vec3_Vector[2].x * vec3_Vector[1].y)
                           / ((vec3_Vector[1].y - vec3_Vector[2].y) * vec3_Vector[0].x + (vec3_Vector[2].x - vec3_Vector[1].x) * vec3_Vector[0].y + vec3_Vector[1].x * vec3_Vector[2].y - vec3_Vector[2].x * vec3_Vector[1].y);
            double beta = ((vec3_Vector[2].y - vec3_Vector[0].y) * i + (vec3_Vector[0].x - vec3_Vector[2].x) * j + vec3_Vector[2].x * vec3_Vector[0].y - vec3_Vector[0].x * vec3_Vector[2].y)
                          / ((vec3_Vector[1].y - vec3_Vector[2].y) * vec3_Vector[0].x + (vec3_Vector[2].x - vec3_Vector[1].x) * vec3_Vector[0].y + vec3_Vector[1].x * vec3_Vector[2].y - vec3_Vector[2].x * vec3_Vector[1].y);
            double gamma = 1 - alpha - beta;
// ------------------------------- //

            // If the point is inside the triangle
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                // Compute the depth (z) of the fragment using barycentric interpolation
                double depth = alpha * vec3_Vector[0].z + beta * vec3_Vector[1].z + gamma * vec3_Vector[2].z;

                // Perform depth test
                if (depth < scene->depthBuffer[i][j]) {
                    // Update the depth buffer
                    scene->depthBuffer[i][j] = depth;

                    // Compute the interpolated color
                    Color c0 = color_Vector[0];
                    Color c1 = color_Vector[1];
                    Color c2 = color_Vector[2];

// ------------------------------- //
                    c0.r = c0.r * alpha;
                    c0.g = c0.g * alpha;
                    c0.b = c0.b * alpha;
                    c1.r = c1.r * beta;
                    c1.g = c1.g * beta;
                    c1.b = c1.b * beta;
                    c2.r = c2.r * gamma;
                    c2.g = c2.g * gamma;
                    c2.b = c2.b * gamma;


                    // Write the interpolated color to the image
                    image[i][j].r = c0.r + c1.r + c2.r;
                    image[i][j].g = c0.g + c1.g + c2.g;
                    image[i][j].b = c0.b + c1.b + c2.b;

// ------------------------------- //
                }
            }
        }
    }
}

bool visible(double den, double num, double &tE, double &tL) {
    if(den > 0){
        double t = num/den;
        if(t > tL) return false;
        else if(t > tE) tE = t;
    }
    else if(den < 0){
        double t = num/den;
        if(t < tE) return false;
        else if(t < tL) tL = t;
    }
    else if(num > 0) return false;
    return true;
}

bool clipping(Vec4 &point1, Vec4 &point2, Color &color1, Color &color2) {
    double tE = 0, tL = 1;
    bool visible_var = false;

// ------------------------------- //
    point1.x /= point1.t;
    point1.y /= point1.t;
    point1.z /= point1.t;
    point1.t /= point1.t;
    
    point2.x /= point2.t;
    point2.y /= point2.t;
    point2.z /= point2.t;
    point2.t /= point2.t;

    double dx = point2.x - point1.x;
    double dy = point2.y - point1.y;
    double dz = point2.z - point1.z;

   

    Color diff_color;
    diff_color.r = color2.r - color1.r;
    diff_color.g = color2.g - color1.g;
    diff_color.b = color2.b - color1.b;

 // ------------------------------- //

    if(visible(dx,-1-point1.x,tE,tL)) {
        if(visible(-dx,point1.x-1,tE,tL)) {
            if(visible(dy,-1-point1.y,tE,tL)) {
                if(visible(-dy,point1.y-1,tE,tL)) {
                    if(visible(dz,-1-point1.z,tE,tL)) {
                        if(visible(-dz,point1.z-1,tE,tL)) {
                            visible_var = true;
                            if(tL < 1) {
                                point2.x = point1.x + tL*dx;
                                point2.y = point1.y + tL*dy;
                                point2.z = point1.z + tL*dz;
                                color2.b = color1.b + diff_color.b * tL;
                                color2.g = color1.g + diff_color.g * tL;
                                color2.r = color1.r + diff_color.r * tL;

                            }
                            if(tE > 0) {
                                point1.x = point1.x + tE*dx;
                                point1.y = point1.y + tE*dy;
                                point1.z = point1.z + tE*dz;
                                color1.b = color1.b + diff_color.b * tE;
                                color1.g = color1.g + diff_color.g * tE;
                                color1.r = color1.r + diff_color.r * tE;
                            }
                        }
                    }
                }
            }
        }
    }
    return visible_var;

}

vector<Vec4> getTransformedPoints(Triangle &triangle, Matrix4 &matrix4, std::vector<Vec3 *> vertices) {
    vector<Vec4> transformedPoints;
    Vec3 first_vec = *vertices[triangle.vertexIds[0] - 1];
    Vec3 second_vec = *vertices[triangle.vertexIds[1]-1];
    Vec3 third_vec = *vertices[triangle.vertexIds[2]-1];


    Vec4 first_vec4 = Vec4(first_vec.x, first_vec.y, first_vec.z, 1, first_vec.colorId);
    Vec4 second_vec4 = Vec4(second_vec.x, second_vec.y, second_vec.z, 1, second_vec.colorId);
    Vec4 third_vec4 = Vec4(third_vec.x, third_vec.y, third_vec.z, 1, third_vec.colorId);

    transformedPoints.push_back(multiplyMatrixWithVec4(matrix4, first_vec4));
    transformedPoints.push_back(multiplyMatrixWithVec4(matrix4, second_vec4));
    transformedPoints.push_back(multiplyMatrixWithVec4(matrix4, third_vec4));

    return transformedPoints;
}

vector<Color> getColorsOfTriangle(Triangle triangle, std::vector<Vec3 *> vertices, Scene* scene) {
    vector<Color> colors;
    Vec3 first_vec = *vertices[triangle.vertexIds[0]-1];
    Vec3 second_vec = *vertices[triangle.vertexIds[1]-1];
    Vec3 third_vec = *vertices[triangle.vertexIds[2]-1];

    colors.push_back(*(scene->colorsOfVertices[first_vec.colorId-1]) );
    colors.push_back(*(scene->colorsOfVertices[second_vec.colorId-1]) );
    colors.push_back(*(scene->colorsOfVertices[third_vec.colorId-1]));
    return colors;
}





void Scene::forwardRenderingPipeline(Camera *camera)
{

    
    
    for (auto mesh : this->meshes) {
        

        // 1 - Modeling Transformation
        Matrix4 modeling_matrix = getTransformationMatrix(mesh, scalings , rotations ,translations );


        // 2 - Camera Transformation
        Matrix4 camera_transformation_matrix = getCameraTransformationMatrix(camera);

        
        // Combine Matrices
        Matrix4 M1 = multiplyMatrixWithMatrix( camera_transformation_matrix , modeling_matrix);


        // 3 - Projection Transformation
        Matrix4  projection_transformation_matrix = getProjectionTransformationMatrix(camera);


        // Combine Matrices
        Matrix4 M2 = multiplyMatrixWithMatrix( projection_transformation_matrix, M1 );


        // 4 - Viewport Transformation
        Matrix4  viewport_transformation_matrix = getViewportTransformationMatrix(camera);



        // 5 - Culling , Clipping , Rasterization
        
        for (auto triangle : mesh->triangles) {
            vector<Vec4> transformedPoints = getTransformedPoints(triangle, M2, this->vertices);
            vector<Color> colors = getColorsOfTriangle(triangle, this->vertices , this);


             Vec3 first_edge = Vec3(transformedPoints[1].x - transformedPoints[0].x, 
                                       transformedPoints[1].y - transformedPoints[0].y, 
                                       transformedPoints[1].z - transformedPoints[0].z, 0);

            Vec3 second_edge = Vec3(transformedPoints[2].x - transformedPoints[0].x, 
                                        transformedPoints[2].y - transformedPoints[0].y, 
                                        transformedPoints[2].z - transformedPoints[0].z, 0);

            Vec3 normal = normalizeVec3(crossProductVec3(first_edge, second_edge));




            if (this->cullingEnabled) {
                if (dotProductVec3(normal, Vec3(transformedPoints[0].x, transformedPoints[0].y, transformedPoints[0].z, 0)) < 0) 
                    continue;
            }
			else{
                if (dotProductVec3(normal, Vec3(0, 0, -1, 0)) < 0) {
                    std::swap(transformedPoints[1], transformedPoints[2]);
                    std::swap(colors[1], colors[2]);
                }
            }


            // 6 - Clipping
            if (mesh->type == 0) {
                Vec4 p1 = transformedPoints[0], p2 = transformedPoints[1], p3 = transformedPoints[2];
                Color c1 = colors[0], c2 = colors[1], c3 = colors[2];
                bool edge1 = clipping(p1, p2, c1, c2);
                bool edge2 = clipping(p2, p3, c2, c3);
                bool edge3 = clipping(p3, p1, c3, c1);
                
                // 7 - Rasterization (Wireframe Mode)
                if (edge1) line_rasterization(p1, p2, c1, c2, viewport_transformation_matrix, camera, image);
                if (edge2) line_rasterization(p2, p3, c2, c3, viewport_transformation_matrix, camera, image);
                if (edge3) line_rasterization(p3, p1, c3, c1, viewport_transformation_matrix, camera, image);
            } 

            // 7 - Rasterization (Solid Mode)
            else
                rasterization(transformedPoints, colors, viewport_transformation_matrix, camera, image, this);
            
        }
    }
}
