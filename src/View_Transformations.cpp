#include "View_Transformations.h"



Matrix4 getCameraTransformationMatrix(Camera *camera){
        Matrix4 rotationMatrix = getIdentityMatrix();
        rotationMatrix.values[0][0] = camera->u.x;
        rotationMatrix.values[0][1] = camera->u.y;
        rotationMatrix.values[0][2] = camera->u.z;
        rotationMatrix.values[1][0] = camera->v.x;
        rotationMatrix.values[1][1] = camera->v.y;
        rotationMatrix.values[1][2] = camera->v.z;
        rotationMatrix.values[2][0] = camera->w.x;
        rotationMatrix.values[2][1] = camera->w.y;
        rotationMatrix.values[2][2] = camera->w.z;

        Matrix4 cameraMatrix = getIdentityMatrix();
        cameraMatrix.values[0][3] = -camera->position.x;
        cameraMatrix.values[1][3] = -camera->position.y;
        cameraMatrix.values[2][3] = -camera->position.z;
        
        Matrix4 camera_transformation_matrix = multiplyMatrixWithMatrix(rotationMatrix, cameraMatrix);

        return camera_transformation_matrix;
}

Matrix4 getProjectionTransformationMatrix(Camera *camera)
{
    Matrix4 projection_transformation_matrix = getIdentityMatrix();

    double dy = (camera->top - camera->bottom); // distance y
    double dx = (camera->right - camera->left); // distance y
    double dz = (camera->far - camera->near); // distance z

    double ny = (camera->top + camera->bottom); // sum y
    double nx = (camera->right + camera->left); // sum y
    double nz = (camera->far + camera->near); // sum z

    if (camera->projectionType == 1) {
        projection_transformation_matrix.values[0][0] = 2 * camera->near / dx;
        projection_transformation_matrix.values[0][2] = (nx/dx);
            
        projection_transformation_matrix.values[1][1] = 2 * camera->near / dy ;
        projection_transformation_matrix.values[1][2] = (ny/dy);

        projection_transformation_matrix.values[2][2] = -(nz/dz);
        projection_transformation_matrix.values[2][3] = -2 * camera->far * camera->near / dz;

        projection_transformation_matrix.values[3][2] = -1;
        projection_transformation_matrix.values[3][3] = 0;
    } 
    else {
        projection_transformation_matrix.values[0][0] = 2/dx;
        projection_transformation_matrix.values[1][1] = 2/dy;
        projection_transformation_matrix.values[2][2] = -(2/dz);
        projection_transformation_matrix.values[0][3] = -(nx/dx);
        projection_transformation_matrix.values[1][3] = -(ny/dy);
        projection_transformation_matrix.values[2][3] = -(nz/dz);
    }

    return projection_transformation_matrix;
}

Matrix4 getViewportTransformationMatrix(Camera *camera){
    Matrix4 viewport_transformation_matrix = getIdentityMatrix();
    viewport_transformation_matrix.values[0][0] = camera->horRes / 2.0;
    viewport_transformation_matrix.values[1][1] = camera->verRes / 2.0;
    viewport_transformation_matrix.values[0][3] = (camera->horRes - 1) / 2.0;
    viewport_transformation_matrix.values[1][3] = (camera->verRes - 1) / 2.0;
    
    return viewport_transformation_matrix;
}


