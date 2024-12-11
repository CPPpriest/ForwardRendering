#include "Rasterize.h"

#include <algorithm>



int makeBetween(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}




int setRange(int val, int valMin, int valMax) {
    if (val < valMin) {
        return valMin;
    } else if (val > valMax) {
        return valMax;
    }
    return val;
}


bool isVisible(double d, double val, double &entr, double &lve) {
    if(d > 0){
        double t = val/d;
        if(t > lve) return false;
        else if(t > entr) entr = t;
    }
    else if(d < 0){
        double t = val/d;
        if(t < entr) return false;
        else if(t < lve) lve = t;
    }
    else if(val > 0) return false;
    return true;
}


bool isAllVisible( double dx, double dy, double dz, Vec4 &sp, double &entr, double &lve)
{
    return isVisible(dx,    -1 - sp.x, entr, lve) &&
           isVisible(-dx,    sp.x - 1, entr, lve) &&
           isVisible(dy,    -1 - sp.y, entr, lve) &&
           isVisible(-dy,    sp.y - 1, entr, lve) &&
           isVisible(dz,    -1 - sp.z, entr, lve) &&
           isVisible(-dz,    sp.z - 1, entr, lve);
}


void interpolatePoint(Vec4 &p, const Vec4 &start, double t, double dx, double dy, double dz) {
    p.x = start.x + t * dx;
    p.y = start.y + t * dy;
    p.z = start.z + t * dz;
}

void interpolateColor(Color &color, const Color &start, const Color &diff, double t) {
    color.r = start.r + diff.r * t;
    color.g = start.g + diff.g * t;
    color.b = start.b + diff.b * t;
}


vector<Vec4> processTrianglePoints(const Triangle &triangle, const Matrix4 &matrix4, const std::vector<Vec3 *> &vertices) {
    vector<Vec4> result;

    for (int vertexId : triangle.vertexIds) {
        Vec3 &vertex = *vertices[vertexId - 1];
        Vec4 vertexVec4(vertex.x, vertex.y, vertex.z, 1, vertex.colorId);
        result.push_back(multiplyMatrixWithVec4(matrix4, vertexVec4));
    }

    return result;
}

vector<Color> getTriangleColors(const Triangle &triangle, const std::vector<Vec3 *> &vertices, const Scene *scene) {
    vector<Color> colors;

    for (int vertexId : triangle.vertexIds) {
        Vec3 &vertex = *vertices[vertexId - 1];
        colors.push_back(*(scene->colorsOfVertices[vertex.colorId - 1]));
    }
    return colors;
}


void drawPixel(std::vector<std::vector<Color>>& image, int x, int y, const Color& cl1, const Color& cl2, double pos, double length) {
    image[x][y].r = makeBetween( (cl1.r * (length - pos) + cl2.r * pos) / length );
    image[x][y].g = makeBetween((cl1.g * (length - pos) + cl2.g * pos) / length);
    image[x][y].b = makeBetween((cl1.b * (length - pos) + cl2.b * pos) / length);
}


void drawSteepNegativeSlope(std::vector<std::vector<Color>>& image, const Vec4& p1, const Vec4& p2, const Color& cl1, const Color& cl2, double dx, double dy) {
    int x = max(p1.x, p2.x);
    int y = min(p1.y, p2.y);

    double valC = 2 * dx - dy;      // current
    double valE = 2 * dx;           // east
    double valNE = 2 * (dx - dy);   // north east

    while (y < max(p1.y, p2.y)) {
        drawPixel(image, x, y, cl1, cl2, abs(y - p2.y), dy);
        if (valC <= 0) {
            valC += valE;
            y++;
        } else {
            valC += valNE;
            x--;
            y++;
        }
    }
}

void drawSteepPositiveSlope(std::vector<std::vector<Color>>& image, const Vec4& p1, const Vec4& p2, const Color& cl1, const Color& cl2, double dx, double dy) {
    int x = min(p1.x, p2.x);
    int y = min(p1.y, p2.y);

    double valC = 2 * dx - dy;
    double valE = 2 * dx;
    double valNE = 2 * (dx - dy);

    while (y < max(p1.y, p2.y)) {
        drawPixel(image, x, y, cl1, cl2, abs(y - p2.y), dy);
        if (valC <= 0) {
            valC += valE;
            y++;
        } else {
            valC += valNE;
            x++;
            y++;
        }
    }
}

void drawShallowNegativeSlope(std::vector<std::vector<Color>>& image, const Vec4& p1, const Vec4& p2, const Color& cl1, const Color& cl2, double dx, double dy) {
    int x = min(p1.x, p2.x);
    int y = max(p1.y, p2.y);

    double valC = 2 * dy - dx;
    double valE = 2 * dy;
    double valNE = 2 * (dy - dx);
    
    while (x < max(p1.x, p2.x)) {
        drawPixel(image, x, y, cl1, cl2, abs(x - p2.x), dx);
        if (valC <= 0) {
            valC += valE;
            x++;
        } else {
            valC += valNE;
            x++;
            y--;
        }
    }
}

void drawShallowPositiveSlope(std::vector<std::vector<Color>>& image, const Vec4& p1, const Vec4& p2, const Color& cl1, const Color& cl2, double dx, double dy) {
    int x = min(p1.x, p2.x);
    int y = min(p1.y, p2.y);
    
    double valC = 2 * dy - dx;
    double valE = 2 * dy; 
    double valNE = 2 * (dy - dx);

    while (x < max(p1.x, p2.x)) {
        drawPixel(image, x, y, cl1, cl2, abs(x - p2.x), dx);
        if (valC <= 0) {
            valC += valE;
            x++;
        } else {
            valC += valNE;
            x++;
            y++;
        }
    }
}


void rasterizeLine( Vec4& p1, Vec4& p2, 
                    Color& cl1, Color& cl2, 
                    Matrix4 matViewport, 
                    std::vector<std::vector<Color>>& image  ) 
{
    double entr = 0;
    double lve = 1;

    // Normalize coordinates
    for (Vec4* point : {&p1, &p2}) {
        point->x /= point->t;
        point->y /= point->t;
        point->z /= point->t;
        point->t = 1;
    }

    // Compute differences
    double d_x = p2.x - p1.x;
    double d_y = p2.y - p1.y;
    double d_z = p2.z - p1.z;
    Color diff = cl2 - cl1;

    // Visibility checks
    if (! isAllVisible( d_x, d_y, d_z, p1, entr , lve) )
        return;

    if (lve < 1) {
        interpolatePoint(p2, p1, lve, d_x, d_y, d_z);
        interpolateColor(cl2, cl1, diff, lve);
    }
    if (entr > 0) {
        interpolatePoint(p1, p1, entr, d_x, d_y, d_z);
        interpolateColor(cl1, cl1, diff, entr);
    }

    Vec4 pv1 = multiplyMatrixWithVec4(matViewport, p1);
    Vec4 pv2 = multiplyMatrixWithVec4(matViewport, p2);

    // Calculate slope and diff
    double slope = (pv2.y - pv1.y) / (pv2.x - pv1.x);
    double dx = abs(pv2.x - pv1.x);
    double dy = abs(pv2.y - pv1.y);

    // start point - end point
    int x = min(pv1.x, pv2.x);
    int y = min(pv1.y, pv2.y);
    int xEnd = max(pv1.x, pv2.x);
    int yEnd = max(pv1.y, pv2.y);

    // Draw line according to slope
    if (slope <= -1) {
        drawSteepNegativeSlope(image, pv1, pv2, cl2, cl1, dx, dy);
        return;
    }
    if (slope >= 1) {
        drawSteepPositiveSlope(image, pv1, pv2, cl2, cl1, dx, dy);
        return;
    }
    if (slope < 0) {
        drawShallowNegativeSlope(image, pv1, pv2, cl2, cl1, dx, dy);
        return;
    }
    if( slope<1 )
    {
        drawShallowPositiveSlope(image, pv1, pv2, cl2, cl1, dx, dy);
        return;
    }
}


void rasterizeSolid(vector<Vec4> points, vector<Color> colour, Matrix4 mat, Camera *camera, std::vector<std::vector<Color>> &image , Scene*scene)
{
    for (int i = 0; i < points.size(); ++i) {
        points[i] = multiplyVec4WithScalar(points[i], 1 / points[i].t);
        points[i] = multiplyMatrixWithVec4(mat, points[i]);
    }

    int xMin = points[0].x;
    int xMax = points[0].x;
    int yMin = points[0].y;
    int yMax = points[0].y;
    
    for (int i = 1; i < points.size(); ++i) 
    {

        if (points[i].x < xMin)
            xMin = points[i].x;
        if (points[i].x > xMax)
            xMax = points[i].x;
        if (points[i].y < yMin) 
            yMin = points[i].y;
        if (points[i].y > yMax)
            yMax = points[i].y;
    }


    xMin = setRange(xMin, 0, camera->horRes - 1);
    yMin = setRange(yMin, 0, camera->verRes - 1);
    xMax = setRange(xMax, 0, camera->horRes - 1);
    yMax = setRange(yMax, 0, camera->verRes - 1);


    vector<Vec3> tp;
    for (int i = 0; i < points.size(); ++i) {
        Vec3 temp = Vec3(points[i].x, points[i].y, points[i].z, points[i].colorId);
        tp.push_back(temp);
    }

    for (int i = xMin; i <= xMax; i++) {
        for (int j = yMin; j <= yMax; j++) {            

            double alpha = (    (tp[1].y - tp[2].y) * i   +   
                                (tp[2].x - tp[1].x) * j   + 
                                (tp[1].x * tp[2].y)       - 
                                (tp[2].x * tp[1].y)
                            )
                           
                           / 
                           
                           (    (tp[1].y - tp[2].y) * tp[0].x    +   
                                (tp[2].x - tp[1].x) * tp[0].y    +   
                                (tp[1].x * tp[2].y)              -   
                                (tp[2].x * tp[1].y)
                            );

            
            double beta =  (    (tp[2].y - tp[0].y) * i   + 
                                (tp[0].x - tp[2].x) * j   + 
                                (tp[2].x * tp[0].y)       - 
                                (tp[0].x * tp[2].y)
                            )

                            / 
                            
                            (
                                (tp[1].y - tp[2].y) * tp[0].x    + 
                                (tp[2].x - tp[1].x) * tp[0].y    + 
                                (tp[1].x * tp[2].y)              - 
                                (tp[2].x * tp[1].y)
                            );
            
            double gamma = 1 - alpha - beta;


            // point inside
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {

                double depth = alpha * tp[0].z + beta * tp[1].z + gamma * tp[2].z;

                // depth test
                if (depth < scene->depthBuffer[i][j]) {
                    // Update depth buffer
                    scene->depthBuffer[i][j] = depth;

                    // interpolated color
                    Color c0 = colour[0];
                    c0 *= alpha;
                    
                    Color c1 = colour[1];
                    c1 *= beta;

                    Color c2 = colour[2];
                    c2 *= gamma;


                    // interpolated color to image
                    image[i][j].r = c0.r + c1.r + c2.r;
                    image[i][j].g = c0.g + c1.g + c2.g;
                    image[i][j].b = c0.b + c1.b + c2.b;

                }
            }
        }
    }
}
