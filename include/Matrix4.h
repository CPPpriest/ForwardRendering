#ifndef __MATRIX4_H__
#define __MATRIX4_H__

#include <initializer_list>


class Matrix4
{
public:
    double values[4][4];

    Matrix4();
    Matrix4(double values[4][4]);
    Matrix4(const Matrix4 &other);
    Matrix4(std::initializer_list<std::initializer_list<double>> init);
    friend std::ostream &operator<<(std::ostream &os, const Matrix4 &m);
};

#endif