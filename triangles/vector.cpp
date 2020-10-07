//
// Created by turfe on 11.04.2020.
//
#include "vector.hpp"
#include <math.h>

Vector3D::Vector3D() : Vector3D(0, 0, 0) {
}

Vector3D::Vector3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {
}

Vector3D::~Vector3D() {
}

double Vector3D::getX() const {
    return x;
}

double Vector3D::getY() const {
    return y;
}

double Vector3D::getZ() const {
    return z;
}

void Vector3D::setX(double x) {
    this->x = x;
}

void Vector3D::setY(double y) {
    this->y = y;
}

void Vector3D::setZ(double z) {
    this->z = z;
}

bool Vector3D::operator==(const Vector3D &v) const {
    double E = 0.000001;
    return (abs(v.x - x) < E) && (abs(v.y - y) < E) && (abs(v.z - z) < E);
}

bool Vector3D::operator!=(const Vector3D &v) const {
    return !(*this == v);
}

Vector3D Vector3D::operator+(const Vector3D &v) const {
    return Vector3D(x + v.x, y + v.y, z + v.z);
}

Vector3D Vector3D::operator-(const Vector3D &v) const {
    return Vector3D(x - v.x, y - v.y, z - v.z);
}

Vector3D Vector3D::operator*(const int a) const {
    return Vector3D(x * a, y * a, z * a);
}

double Vector3D::operator*(const Vector3D &v) const {
    return x * v.x + y * v.y + z * v.z;
}


Vector3D operator*(int a, const Vector3D &v) {
    return Vector3D(v.getX() * a, v.getY() * a, v.getZ() * a);
}

std::ostream &operator<<(std::ostream &os, const Vector3D &v) {
    os << "(" << v.getX() << "; " << v.getY() << "; " << v.getZ() << ")";
    return os;
}

std::istream &operator>>(std::istream &is, Vector3D &v) {
    double x, y, z;
    is >> x >> y >> z;
    v.setX(x);
    v.setY(y);
    v.setZ(z);
    return is;
}

Vector3D::Vector3D( const triangle_geometry::Point& point1, const triangle_geometry::Point& point2):
x(point1.get_x() - point2.get_x()), y(point1.get_y() - point2.get_y()), z(point1.get_z() - point2.get_z()){}