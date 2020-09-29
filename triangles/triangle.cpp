#include "triangle.h"

//it is the maaaaaain function of this problem
bool triangle_geometry::triangle_handler::is_intersect(const triangle &tr1, const triangle &tr2) const {
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
    Plane plane_1 = Plane( tr1);
    int sign1, sign2, sign3;
    sign1 = triangle_handler::sign_of_dist(plane_1, tr2.getVertice(1));
    sign2 = triangle_handler::sign_of_dist(plane_1, tr2.getVertice(2));
    sign3 = triangle_handler::sign_of_dist(plane_1, tr2.getVertice(3));
    if ( sign1 * sign2 > 0 && sign1 * sign3 > 0) return false;//треугольник полностью с одной стороны от плоскости
    //else continue
    Plane plane_2 = Plane( tr2);
    if ( triangle_handler::is_coincident( plane_1, plane_2, e)) {
        if (plane_1.get_d() == plane_2.get_d()) {
            /*
             * пересечение треугольников, лежащих в одной плоскости
             */
        }return false;//треугольники в параллельных плоскостях
    }
    sign1 = triangle_handler::sign_of_dist(plane_2, tr1.getVertice(1));
    sign2 = triangle_handler::sign_of_dist(plane_2, tr1.getVertice(2));
    sign3 = triangle_handler::sign_of_dist(plane_2, tr1.getVertice(3));
    if ( sign1 * sign2 > 0 && sign1 * sign3 > 0) return false;//треугольник полностью с одной стороны от плоскости
    /*
     * вычисление линии пересечения плоскостей
     */
    /*
     * вычисление "проекций" на линию пересечения
     */
    /*
     * если нет пересечений, то треуглоьники не пересекаются,
     * иначе - пересекаются, и не надо дальше считать
     */
}

int triangle_geometry::triangle_handler::sign_of_dist(const Plane &plane, const Point &point) const {
    double dist = triangle_handler::sc_pr(plane.get_norm_vec(), point) + plane.get_d();
    if ( dist < 0) return -1;
    if ( dist > 0) return 1;
    else return 0;
}

triangle_geometry::Plane::Plane() = default;

triangle_geometry::Plane::Plane(const Point& norm_vec, const Point &point){
    _norm_vec = norm_vec;
    _point = point;
    _d = triangle_handler::sc_pr( norm_vec, point);
}

triangle_geometry::Plane::Plane(const triangle_geometry::triangle &tr): Plane( triangle_handler::vec_mul(tr.getVertice(1) - tr.getVertice(2),
                                                                          tr.getVertice(2) - tr.getVertice(3)), tr.getVertice(3) ) {}

triangle_geometry::Point triangle_geometry::Plane::get_norm_vec() const{
    return _norm_vec;
}

double triangle_geometry::Plane::get_d() const{
    return _d;
}

triangle_geometry::triangle::triangle() = default;

bool triangle_geometry::triangle::is_degenerate() const{//true if degenerate
    Point norm_vec = triangle_handler::vec_mul( _p1 - _p2, _p2 - _p3 );
    if ( norm_vec.get_x() == 0 && norm_vec.get_y() == 0 && norm_vec.get_z() == 0) return true;
    return false;
}

triangle_geometry::Point::Point(const double &x, const double &y, const double &z): _x(x), _y(y), _z(z) {}

triangle_geometry::Point::Point():Point(0.0, 0.0, 0.0)  {}

triangle_geometry::Point::Point(const Point &source) = default;

triangle_geometry::Point triangle_geometry::Point::operator-(const Point &source) const {
    return Point(_x - source._x, _y - source._y, _y - source._z);
}

double triangle_geometry::Point::get_x() const {
    return _x;
}

double triangle_geometry::Point::get_y() const {
    return _y;
}

double triangle_geometry::Point::get_z() const {
    return _z;
}

triangle_geometry::Point triangle_geometry::triangle::getVertice(const int &i) const {
    if ( i == 1) return _p1;
    if ( i == 2) return _p2;
    return _p3;
}

triangle_geometry::triangle::triangle(const triangle_geometry::Point &p1, const triangle_geometry::Point &p2,
                                      const triangle_geometry::Point &p3): _p1(p1), _p2(p2), _p3(p3) {}

triangle_geometry::Point triangle_geometry::triangle_handler::vec_mul(const Point &p1, const Point &p2){
    return Point(p1.get_y() * p2.get_z() - p1.get_z() * p2.get_y(), p1.get_z() * p2.get_x() - p1.get_x() * p2.get_z(),
                 p1.get_x() * p2.get_y() - p1.get_y() * p2.get_x());
}

double triangle_geometry::triangle_handler::sc_pr(const Point &p1, const Point &p2){
    return ( p1.get_x() * p2.get_x() + p1.get_y() * p2.get_y() + p1.get_z() * p2.get_z());
}

bool triangle_geometry::triangle_handler::is_coincident(const triangle_geometry::Plane &p1,
                                                        const triangle_geometry::Plane &p2, const double &e) const {
    if ( (p1.get_norm_vec().get_x() - e) <= p2.get_norm_vec().get_x() && (p1.get_norm_vec().get_x() + e) >= p2.get_norm_vec().get_x()
    && (p1.get_norm_vec().get_y() - e) <= p2.get_norm_vec().get_y() && (p1.get_norm_vec().get_y() + e) >= p2.get_norm_vec().get_y()
    && (p1.get_norm_vec().get_z() - e) <= p2.get_norm_vec().get_z() && (p1.get_norm_vec().get_z() + e) >= p2.get_norm_vec().get_z()){
        return true;
    }
    return false;
}
