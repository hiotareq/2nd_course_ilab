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
        Point(const Point& source);
        Point operator-( const Point& p2) const;
        double get_x() const;
        double get_y() const;
        double get_z() const;
        double operator*( const Point& source) const;//скалярное умножение, точки == векторы
        Point operator%( const Point& source) const;//векторное умножение, точки == векторы
        Point operator+(const Point& source) const;//точка == вектор
        Point operator/(const double& source) const;
        Point operator*( const int& source) const;//точка == вектор
        Point operator*( const double& source) const;//точка == вектор
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
        Point _direction;//направляющий вектор прямой
        Point _point;//точка на прямой
    public:
        Point GetPoint();
        Point Getdirection();
        Line(const Point& direction, const Point& point);
        bool IsOverlap( const Point& p1,const Point& p2,const Point& p3,const Point& p4,const Point& p5,const Point& p6);
    };

    class Plane{
    private:
        Point _norm_vec, _point;
        double _d;
    public:
        double get_d() const;
        Point get_norm_vec() const;
        Plane();
        Plane( const Point& norm_vec, const Point &point);
        Plane( const triangle_geometry::triangle & tr);
    };

    class triangle_handler{
    private:
        double e = 1e-3;
    public:
        bool is_intersect( const triangle& tr1,const triangle& tr2) const;
        int sign_of_dist( const Plane &plane, const Point& point ) const;//знак расстояния от точки до плоскости
        bool is_coincident( const Plane& p1, const Plane& p2) const;//true, если совпадают
        bool is_intersect2D( const triangle& tr1, const triangle& tr2) const;//true, если пересекаются
        int GetMiddleIndex( const int& i0, const int& i1) const;
        int GetExtremeIndex( const triangle& tr, const Point& point) const;
        Point GetNormVector( const triangle& tr, const Point& p1, const Point& p2) const;
        Line GetLine(const Plane& p1, const Plane& p2) const;
    };
}

#endif //INC_2ND_TERM_ILAB_TRIANGLE_H
