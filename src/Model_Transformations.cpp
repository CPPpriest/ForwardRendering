#include "Model_Transformations.h"




Matrix4 getModelMat(Mesh *mesh, std::vector<Scaling *>  scalings , std::vector<Rotation *> rotations , std::vector<Translation *> translations ) {
    
    Matrix4 resMat = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; i++) {
        int id = mesh->transformationIds[i];    
        --id;

        switch (mesh->transformationTypes[i]) {
            case 's': 
            {
				double sx = scalings[id]->sx;
				double sy = scalings[id]->sy;
				double sz = scalings[id]->sz;

				Matrix4 matScale = getIdentityMatrix();
				matScale.values[0][0] = sx;
				matScale.values[1][1] = sy;
				matScale.values[2][2] = sz;
				
                resMat = multiplyMatrixWithMatrix(matScale, resMat);		
                break;
            }
            case 'r': 
            {
                double alpha = rotations[id]->angle;
                double ux = rotations[id]->ux;
                double uy = rotations[id]->uy;
                double uz = rotations[id]->uz;

                Matrix4 matRotate = getIdentityMatrix();

                double c = cos(alpha * M_PI / 180.0);
                double s = sin(alpha * M_PI / 180.0);
                double t = 1 - c;

                double values[3][3] = {
                    {t * ux * ux + c,        t * ux * uy - s * uz, t * ux * uz + s * uy},
                    {t * ux * uy + s * uz,   t * uy * uy + c,      t * uy * uz - s * ux},
                    {t * ux * uz - s * uy,   t * uy * uz + s * ux, t * uz * uz + c}
                };

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        matRotate.values[i][j] = values[i][j];
                    }
                }

                resMat = multiplyMatrixWithMatrix(matRotate, resMat);
                break;
            }
            case 't': 
            {
                double tx = translations[id]->tx;
                double ty = translations[id]->ty;
                double tz = translations[id]->tz;

                Matrix4 matTranslate = getIdentityMatrix();

                matTranslate.values[0][3] = tx;
                matTranslate.values[1][3] = ty;
                matTranslate.values[2][3] = tz;

                resMat = multiplyMatrixWithMatrix(matTranslate, resMat);
                break;
            }
        }
    }

    return resMat;
}



