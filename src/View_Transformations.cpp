#include "View_Transformations.h"



Matrix4 getCameraMat(Camera *camera){
        
        Matrix4 matRotate = getIdentityMatrix();

        matRotate.values[0][0] = camera->u.x;
        matRotate.values[0][1] = camera->u.y;
        matRotate.values[0][2] = camera->u.z;
        matRotate.values[1][0] = camera->v.x;
        matRotate.values[1][1] = camera->v.y;
        matRotate.values[1][2] = camera->v.z;
        matRotate.values[2][0] = camera->w.x;
        matRotate.values[2][1] = camera->w.y;
        matRotate.values[2][2] = camera->w.z;

        Matrix4 matCamera = getIdentityMatrix();
        matCamera.values[0][3] = -camera->position.x;
        matCamera.values[1][3] = -camera->position.y;
        matCamera.values[2][3] = -camera->position.z;
        
        return multiplyMatrixWithMatrix(matRotate, matCamera);
}

Matrix4 getProjectionMat(Camera *camera)
{
    Matrix4 matProjection = getIdentityMatrix();

    double dy = (camera->top - camera->bottom); // distance y
    double dx = (camera->right - camera->left); // distance y
    double dz = (camera->far - camera->near); // distance z

    double ny = (camera->top + camera->bottom); // sum y
    double nx = (camera->right + camera->left); // sum y
    double nz = (camera->far + camera->near); // sum z

    // PERSPECTIVE_PROJECTION
    if (camera->projectionType == PERSPECTIVE_PROJECTION) {

        matProjection.values[0][0] = 2 * camera->near / dx;
        matProjection.values[0][2] = (nx/dx);
            
        matProjection.values[1][1] = 2 * camera->near / dy ;
        matProjection.values[1][2] = (ny/dy);

        matProjection.values[2][2] = -(nz/dz);
        matProjection.values[2][3] = -2 * camera->far * camera->near / dz;

        matProjection.values[3][2] = -1;
        matProjection.values[3][3] = 0;
    } 
    // ORTOGRAPHIC_PROJECTION
    else {
        matProjection.values[0][0] = 2/dx;
        matProjection.values[1][1] = 2/dy;
        matProjection.values[2][2] = -(2/dz);

        matProjection.values[0][3] = -(nx/dx);
        matProjection.values[1][3] = -(ny/dy);
        matProjection.values[2][3] = -(nz/dz);
    }

    return matProjection;
}

Matrix4 getViewportMat(Camera *camera){
    Matrix4 matViewport = getIdentityMatrix();
    matViewport.values[0][0] = camera->horRes / 2.0;
    matViewport.values[1][1] = camera->verRes / 2.0;
    matViewport.values[0][3] = (camera->horRes - 1) / 2.0;
    matViewport.values[1][3] = (camera->verRes - 1) / 2.0;
    
    return matViewport;
}


