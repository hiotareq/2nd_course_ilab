#ifndef INC_2ND_TERM_ILAB_TRIANGLE_H
#define INC_2ND_TERM_ILAB_TRIANGLE_H

#include <iostream>
#include <vector>

namespace triangle_geometry{
    struct Point{
        double _x,_y,_z;
        Point(const double &x, &y, &z);
        Point();
        Point(const Point& source)
        Point operator-( const Point& p2) const;
    };

    struct Plane{
    private:
        Point _norm_vec, _point;
        double _d;
    public:
        double get_d();
        Point get_norm_vec();
        Plane();
        Plane(const Point& norm_vec, &point);
    };

    class triangle{
    private:
        std::vector<Point>[3] vertices;
    public:
        triangle();
        triangle(const Point &p1, &p2, &p3);
        bool is_degenerate() const;//true if degenerate
        Point get_vertice(const int &i) const;
    };

    class triangle_handler{
    private:
    public:
        bool is_intersect(triangle& tr1, triangle& tr2) const;
        bool vec_mul(const Point &p1, &p2 ) const;//true if !=0, передаются точки как векторы
        Point vec_mul(const Point &p1, &p2 ) const;//снова точки как векторы
        Plane cal_plane(const triangle& tr) const;
        int sign_of_dist( const Plane &plane, const Point& point ) const;
        double sc_pr(const Point& p1, &p2) const;
    };
}

#endif //INC_2ND_TERM_ILAB_TRIANGLE_H
