//
// Created by turfe on 11.04.2020.
//

#include <iostream>
#include "triangle.hpp"

#ifndef VECTOR_H
#define VECTOR_H

class Vector3D {
private:
    double x, y, z;
public:
    Vector3D();

    Vector3D(double _x, double _y, double _z);

    Vector3D(const triangle_geometry::Point& point);

    Vector3D(const triangle_geometry::Point& point1, const triangle_geometry::Point& point2);

    ~Vector3D();

    double getX() const;

    double getY() const;

    double getZ() const;

    void setX(double x);

    void setY(double y);

    void setZ(double z);

    bool operator==(const Vector3D &v) const;

    bool operator!=(const Vector3D &v) const;

    Vector3D operator+(const Vector3D &v) const;

    Vector3D operator-(const Vector3D &v) const;

    Vector3D operator*(const int a) const;

    double operator*(const Vector3D &v) const;

    double length() const;

};

Vector3D operator*(int a, const Vector3D &v);

std::ostream& operator<<(std::ostream& os, const Vector3D &v);

std::istream& operator>>(std::istream& is, Vector3D &v);

#endif //LINAL_VECTOR_H
