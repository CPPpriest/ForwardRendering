#include "View_Transformations.h"

Matrix4 getCameraMat(Camera *camera) {
    // rotation matrix
    Matrix4 matRotate = {
        {camera->u.x, camera->u.y, camera->u.z, 0},
        {camera->v.x, camera->v.y, camera->v.z, 0},
        {camera->w.x, camera->w.y, camera->w.z, 0},
        {0, 0, 0, 1}
    };

    // translation matrix
    Matrix4 matTranslate = getIdentityMatrix();
    matTranslate.values[0][3] = -camera->position.x;
    matTranslate.values[1][3] = -camera->position.y;
    matTranslate.values[2][3] = -camera->position.z;

    // Combine matrices
    return multiplyMatrixWithMatrix(matRotate, matTranslate);
}

Matrix4 getProjectionMat(Camera *camera) {
    Matrix4 matProjection = getIdentityMatrix();

    double dx = camera->right - camera->left;
    double dy = camera->top - camera->bottom;
    double dz = camera->far - camera->near;

    double nx = camera->right + camera->left;
    double ny = camera->top + camera->bottom;
    double nz = camera->far + camera->near;

    if (camera->projectionType == PERSPECTIVE_PROJECTION) {
        // Perspective
        matProjection.values[0][0] = 2 * camera->near / dx;
        matProjection.values[0][2] = nx / dx;

        matProjection.values[1][1] = 2 * camera->near / dy;
        matProjection.values[1][2] = ny / dy;

        matProjection.values[2][2] = -nz / dz;
        matProjection.values[2][3] = -2 * camera->far * camera->near / dz;

        matProjection.values[3][2] = -1;
        matProjection.values[3][3] = 0;
    } else {
        // Orthographic
        matProjection.values[0][0] = 2 / dx;
        matProjection.values[1][1] = 2 / dy;
        matProjection.values[2][2] = -2 / dz;

        matProjection.values[0][3] = -nx / dx;
        matProjection.values[1][3] = -ny / dy;
        matProjection.values[2][3] = -nz / dz;
    }

    return matProjection;
}

Matrix4 getViewportMat(Camera *camera) {
    return {
        {camera->horRes / 2.0, 0, 0, (camera->horRes - 1) / 2.0},
        {0, camera->verRes / 2.0, 0, (camera->verRes - 1) / 2.0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
}
