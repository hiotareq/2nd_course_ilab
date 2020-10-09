#ifndef INC_2ND_TERM_ILAB_TRIANGLE_H
#define INC_2ND_TERM_ILAB_TRIANGLE_H

#include <iostream>
#include <vector>
#include <map>

namespace triangle_geometry{
    class Point{
        double _x,_y,_z;
    public:
        Point( const double &x,const double &y,const double &z);
        Point();
        Point(const Vector3D& vector);
        Point(const Point& p) = default;
        Point& operator=(const Point& p) = default;
        double get_x() const;
        double get_y() const;
        double get_z() const;
        double length() const;
    };

    class triangle{
    private:
        Point _p1, _p2, _p3;
    public:
        friend Point;
        triangle();
        triangle( const Point &p1, const Point &p2, const Point &p3);
        bool is_degenerate() const;//true if degenerate
        Point getVertice( const int &i) const;
    };

    class Line{
    private:
        Vector3D _direction;//направляющий вектор прямой
        Point _point;//точка на прямой
    public:
        Point GetPoint() const;
        Vector3D Getdirection() const;
        Line(const Vector3D& direction, const Point& point);
        bool IsOverlap( const Point& p1,const Point& p2,const Point& p3,const Point& p4,const Point& p5,const Point& p6);
    };

    class Plane{
    private:
        Point _point;
        Vector3D normal;
        double _d;
    public:
        double get_d() const;
        Vector3D get_normal() const;
        Plane();
        Plane( const Vector3D& norm_vec, const Point &point);
        Plane( const triangle_geometry::triangle & tr);
    };

    class triangle_handler{
    private:
        double e = 1e-3;
    public:
        bool is_intersect( const triangle& tr1,const triangle& tr2) const;
        int sign_of_dist( const Plane &plane, const Point& point ) const;//знак расстояния от точки до плоскости
        bool is_coincident( const Plane& p1, const Plane& p2, const Point& point1, const Point& point2) const;//true, если совпадают
        bool is_intersect2D( const triangle& tr1, const triangle& tr2) const;//true, если пересекаются
        int GetMiddleIndex( const int& i0, const int& i1) const;
        int GetExtremeIndex( const triangle& tr, const Vector3D& point) const;
        Vector3D GetNormVector( const triangle& tr, const Point& p1, const Point& p2) const;
        Line GetLine(const Plane& p1, const Plane& p2) const;
        Point IntersectionEdgeLine(const Point& PointFromTriangle,const Line& line, const Vector3D& EdgeDir) const;
    };

    class Vector3D {
    private:
        double _x, _y, _z;
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

        Vector3D operator%(const Vector3D &v) const;

    };

Point operator+(const Vector3D& v, const Point& p);

Vector3D operator*(const int& a, const Vector3D &v);

std::ostream& operator<<(std::ostream& os, const Vector3D &v);

std::istream& operator>>(std::istream& is, Vector3D &v);
}

#endif //INC_2ND_TERM_ILAB_TRIANGLE_H
