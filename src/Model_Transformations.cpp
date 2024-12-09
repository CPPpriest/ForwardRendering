#include "Model_Transformations.h"




Matrix4 getTransformationMatrix(Mesh *mesh, std::vector<Scaling *>  scalings , std::vector<Rotation *> rotations , std::vector<Translation *> translations ) {
    Matrix4 transformation_matrix = getIdentityMatrix();
    for (int i = 0; i < mesh->numberOfTransformations; i++) {
        int transformationId = mesh->transformationIds[i];
        char transformationType = mesh->transformationTypes[i];
        switch (transformationType) {
            case 's': {
				// Assuming you have access to sx, sy, sz for the current scaling
				double sx = scalings[transformationId - 1]->sx;
				double sy = scalings[transformationId - 1]->sy;
				double sz = scalings[transformationId - 1]->sz;

				// Create an identity matrix for scaling
				Matrix4 scalingMatrix = getIdentityMatrix();
				scalingMatrix.values[0][0] = sx;
				scalingMatrix.values[1][1] = sy;
				scalingMatrix.values[2][2] = sz;

				// Multiply the transformation_matrix with scalingMatrix and store the result in transformation_matrix
				transformation_matrix = multiplyMatrixWithMatrix(scalingMatrix, transformation_matrix);
				break;
                // this->scalings[transformationId - 1]->applyScaling(transformation_matrix);
                // break;
                }
            case 'r': {
				 // Assuming you have access to the angle and axis components (ux, uy, uz) for the current rotation
                double angle = rotations[transformationId - 1]->angle;
                double ux = rotations[transformationId - 1]->ux;
                double uy = rotations[transformationId - 1]->uy;
                double uz = rotations[transformationId - 1]->uz;

                // Create an identity matrix for rotation
                Matrix4 rotationMatrix = getIdentityMatrix();

                // Calculate cosine, sine, and (1 - cosine)
                double c = cos(angle * M_PI / 180.0);
                double s = sin(angle * M_PI / 180.0);
                double t = 1 - c;

                // Fill the rotation matrix
                rotationMatrix.values[0][0] = t * ux * ux + c;
                rotationMatrix.values[0][1] = t * ux * uy - s * uz;
                rotationMatrix.values[0][2] = t * ux * uz + s * uy;
                rotationMatrix.values[1][0] = t * ux * uy + s * uz;
                rotationMatrix.values[1][1] = t * uy * uy + c;
                rotationMatrix.values[1][2] = t * uy * uz - s * ux;
                rotationMatrix.values[2][0] = t * ux * uz - s * uy;
                rotationMatrix.values[2][1] = t * uy * uz + s * ux;
                rotationMatrix.values[2][2] = t * uz * uz + c;

                // Multiply the transformation_matrix with rotationMatrix and store the result in transformation_matrix
                transformation_matrix = multiplyMatrixWithMatrix(rotationMatrix, transformation_matrix);
                break;
                // this->rotations[transformationId - 1]->applyRotation(transformation_matrix);
                // break;
                }
            case 't': {
                 // Assuming you have access to the translation values (tx, ty, tz) for the current transformation
                double tx = translations[transformationId - 1]->tx;
                double ty = translations[transformationId - 1]->ty;
                double tz = translations[transformationId - 1]->tz;

                // Create an identity matrix for translation
                Matrix4 translationMatrix = getIdentityMatrix();

                // Set the translation values
                translationMatrix.values[0][3] = tx;
                translationMatrix.values[1][3] = ty;
                translationMatrix.values[2][3] = tz;

                // Multiply the transformation_matrix with translationMatrix and store the result in transformation_matrix
                transformation_matrix = multiplyMatrixWithMatrix(translationMatrix, transformation_matrix);
                break;
                // this->ranslations[transformationId - 1]->applyTranslation(transformation_matrix);
                // break;
                }
            default:
                break;
        }
    }

    return transformation_matrix;
}



