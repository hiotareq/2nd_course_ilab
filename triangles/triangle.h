#ifndef INC_2ND_TERM_ILAB_TRIANGLE_H
#define INC_2ND_TERM_ILAB_TRIANGLE_H

#include <iostream>
#include <vector>

namespace triangle_geometry{
    class Point{
        double _x,_y,_z;
    public:
        Point(const double &x,const double &y,const double &z);
        Point();
        Point(const Point& source);
        Point operator-( const Point& p2) const;
        double get_x() const;
        double get_y() const;
        double get_z() const;
    };

    class triangle{
    private:
        Point _p1, _p2, _p3;
    public:
        triangle();
        triangle(const Point &p1, const Point &p2, const Point &p3);
        bool is_degenerate() const;//true if degenerate
        Point getVertice(const int &i) const;
    };

    class Plane{
    private:
        Point _norm_vec, _point;
        double _d;
    public:
        double get_d() const;
        Point get_norm_vec() const;
        Plane();
        Plane(const Point& norm_vec, const Point &point);
        Plane( const triangle_geometry::triangle & tr);
    };

    class triangle_handler{
    private:
        double e = 0.00001;
    public:
        bool is_intersect(const triangle& tr1,const triangle& tr2) const;
        static Point vec_mul(const Point &p1, const Point &p2 );//точки = векторы
        int sign_of_dist( const Plane &plane, const Point& point ) const;
        static double sc_pr(const Point& p1, const Point& p2);//точки = векторы
        bool is_coincident( const Plane& p1, const Plane& p2, const double& e) const;//true, если совпадают
    };
}

#endif //INC_2ND_TERM_ILAB_TRIANGLE_H
