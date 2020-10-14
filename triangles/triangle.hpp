#ifndef INC_2ND_TERM_ILAB_TRIANGLE_HPP
#define INC_2ND_TERM_ILAB_TRIANGLE_HPP

#include <iostream>
#include <vector>
#include <map>

namespace triangle_geometry{
    class Point{
    float _x,_y,_z;
    public:
        Point( const float &x,const float &y,const float &z);
        Point();

        Point(const Point& p) = default;
        Point& operator=(const Point& p) = default;

        float get_x() const;
        float get_y() const;
        float get_z() const;

        bool operator==(const Point& p) const;

        void set_x(const float& newx);
        void set_y(const float& newy);
        void set_z(const float& newz);

        void set_point(const float& newx, const float& newy, const float& newz);
    };

    class Vector3D {
    private:
        float _x, _y, _z;
    public:
        Vector3D();

        Vector3D(float _x, float _y, float _z);

        Vector3D(const triangle_geometry::Point& point);

        Vector3D(const triangle_geometry::Point& point1, const triangle_geometry::Point& point2);

        ~Vector3D();

        float getX() const;

        float getY() const;

        float getZ() const;

        void setX(float x);

        void setY(float y);

        void setZ(float z);

        bool operator==(const Vector3D &v) const;

        bool operator!=(const Vector3D &v) const;

        Vector3D operator+(const Vector3D &v) const;

        Vector3D operator-(const Vector3D &v) const;

        Vector3D operator*(const int a) const;

        Vector3D operator*(const float a) const;

        float operator*(const Vector3D &v) const;

        float length() const;

        Vector3D operator%(const Vector3D &v) const;

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

        double get_max_x() const;
        double get_max_y() const;
        double get_max_z() const;
        double get_min_x() const;
        double get_min_y() const;
        double get_min_z() const;

        void set_triangle(const float& x1, const float& y1, const float& z1, const float& x2, const float& y2, const float& z2, const float& x3, const float& y3, const float& z3);

        float get_max_coord();
        float get_min_coord();

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
        float _d;
    public:
        float get_d() const;
        Vector3D get_normal() const;
        Plane();
        Plane( const Vector3D& norm_vec, const Point &point);
        Plane( const triangle_geometry::triangle & tr);
    };

bool is_intersect( const triangle& tr1,const triangle& tr2);

int sign_of_dist( const Plane &plane, const Point& point );//знак расстояния от точки до плоскости

bool is_coincident( const Plane& p1, const Plane& p2, const Point& point1, const Point& point2);//true, если совпадают

bool is_intersect2D( const triangle& tr1, const triangle& tr2);//true, если пересекаются

int GetMiddleIndex( const int& i0, const int& i1);

int GetExtremeIndex( const triangle& tr, const Vector3D& point);

Vector3D GetNormVector( const triangle& tr, const Point& p1, const Point& p2);

Line GetLine(const Plane& p1, const Plane& p2);

Point IntersectionEdgeLine(const Point& PointFromTriangle,const Line& line, const Vector3D& EdgeDir);

Point MakePointFromVector(const Vector3D& v);

Vector3D MakeVectorFromPoint(const Point& p);

std::vector<triangle_geometry::Point> DefinePoints( const Point& p1, const Point& p2, const Point& p3, const Point& p4);

std::vector<triangle_geometry::Point> MakeVector(const Point& most_far, const Point& p1, const Point& p2, const Point& p3);

Point MostFarPoint(const Point& point_to_cmp, const Point& p1, const Point& p2, const Point& p3 );

Point operator+(const Point& p, const Vector3D& v);

std::ostream& operator<<(std::ostream& os, const Vector3D &v);

std::istream& operator>>(std::istream& is, Vector3D &v);

std::istream& operator>>(std::istream& is, triangle& tr);
}

triangle_geometry::Vector3D operator*(const int& a, const triangle_geometry::Vector3D& v);

triangle_geometry::Vector3D operator*(const float& a, const triangle_geometry::Vector3D& v);

#endif //INC_2ND_TERM_ILAB_TRIANGLE_H
