#include <iostream>
#include "triangle.h"
#include <cmath>

//it is the maaaaaain function of this problem
bool triangle_geometry::triangle_handler::is_intersect(const triangle &tr1, triangle &tr2) const {
    if ( tr1.is_degenerate() ){
        /*
         * удаление вырожденного треугольника из базы данных
         */
    }
    if ( tr2.is_degenerate()){
        /*
         * удаление вырожденного треугольника из базы данных
         */
    }
    Plane plane_1 = cal_plane(tr1 );
    int sign1, sign2, sign3;
    sign1 = triangle_handler::sign_of_dist(plane_1, tr2.get_vertice(1));
    sign2 = triangle_handler::sign_of_dist(plane_1, tr2.get_vertice(2));
    sign3 = triangle_handler::sign_of_dist(plane_1, tr2.get_vertice(3));
    if ( sign1 * sign2 > 0 && sign1 * sign3 > 0) return false;//одинаковые знаки у всех расстояний
    //else continue
    Plane plane_2 = cal_plane(tr2 );
    if ( plane_1.get_norm_vec()._x/plane_2.get_norm_vec().x == plane_1.get_norm_vec()._y/plane_2.get_norm_vec()._y
        && plane_1.get_norm_vec()._y/plane_2.get_norm_vec()._y == plane_1.get_norm_vec()._z/plane_2.get_norm_vec()._z) {
        if (plane_1.get_d() == plane_2.get_d()) {
            /*
             * пересечение треугольников, лежащих в одной плоскости
             */
        } else return false;//треугольники в параллельных плоскостях
    }

}

int triangle_geometry::triangle_handler::sign_of_dist(const Plane &plane, const Point &point) const {
    double dist = triangle_handler::sc_pr(plane.get_norm_vec(), point) + plane.get_d();
    if ( dist < 0) return -1;
    if ( dist > 0) return 1;
    else return 0;
}

Plane triangle_geometry::triangle_handler::cal_plane(const triangle &tr) const {
    return Plane(vec_mul(tr.get_vertice(1) - tr.get_vertice(2), tr.get_vertice(2) - tr.get_vertice(3)), tr.get_vertice(3));
}

triangle_geometry::Plane::Plane():Plane(0, 0) {}

triangle_geometry::Plane::Plane(const Point& norm_vec, &point):_norm_vec(norm_vec), _point(point),
                                                                _d(triangle_handler::sc_pr(norm_vec, point)) {}

Point triangle_geometry::Plane::get_norm_vec() {
    return _norm_vec;
}

double triangle_geometry::Plane::get_d() {
    return d;
}

triangle_geometry::triangle::triangle(const Point &p1, int &p2, int &p3):(vertices[0])(p1),
                                                                (vertices[1])(p2), (vertices[2])(p3) {}

triangle_geometry::triangle::triangle():triangle(0.0, 0.0, 0.0){}

bool triangle_geometry::triangle::is_degenerate() const{//true if degenerate
    Point norm_vec = triangle_handler::vec_mul(vertices[0] - vertices[1], vertices[1] - vertices[2] );
    if ( norm_vec._x == 0 && norm_vec._y == 0 && norm_vec._z == 0) return true;
    return false;
}

triangle_geometry::Point::Point(const double &x, int &y, int &z):_x(x), _y(y), _z(z) {}

triangle_geometry::Point::Point():Point(0.0, 0.0, 0.0)  {}

triangle_geometry::Point::Point(const Point &source):_x(source._x), _y(source._y), _z(source._z) {}

Point triangle_geometry::Point::operator-(const Point &p2) const {
    return Point(x - p2._x, y - p2._y, z - p2._z);
}

Point triangle_geometry::triangle::get_vertice(const int &i) const {
    return vertices[i-1];
}

Point triangle_geometry::triangle_handler::vec_mul(const Point &p1, &p2) const {
    return Point(p1._y * p2._z - p1._z * p2._y, p1._z * p2._x - p1._x * p2._z, p1._x * p2._y - p1._y * p2._x);
}

double triangle_geometry::triangle_handler::sc_pr(const Point &p1, int &p2) const {
    return ( p1._x * p2._x + p1._y * p2._y + p1._z * p2._z);
}