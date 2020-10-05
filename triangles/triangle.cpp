#include "triangle.h"
#include <cmath>


//TO-DO: точка пересечения стороны треугольника и прямой пересечения плоскостей, пересечение отрезков, 
//этими прямыми

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
    if ( triangle_handler::is_coincident( plane_1, plane_2, tr1.getVertice(1), tr2.getVertice(1) )) {
        if (plane_1.get_d() == plane_2.get_d()) {
            if ( is_intersect2D(tr1, tr2)) return true;
            return false;
        }return false;//треугольники в параллельных плоскостях
    }
    sign1 = triangle_handler::sign_of_dist(plane_2, tr1.getVertice(1));
    sign2 = triangle_handler::sign_of_dist(plane_2, tr1.getVertice(2));
    sign3 = triangle_handler::sign_of_dist(plane_2, tr1.getVertice(3));
    if ( sign1 * sign2 > 0 && sign1 * sign3 > 0) return false;//треугольники полностью с одной стороны от плоскости
    Line IntersectLine = triangle_handler::GetLine(plane_1, plane_2);
    Point p11, p12, p13, p21, p22, p23;
    p11 = IntersectLine.GetPoint() + IntersectLine.Getdirection() * ( (tr1.getVertice(0) - IntersectLine.GetPoint()) * IntersectLine.Getdirection()) /
            (IntersectLine.Getdirection() * IntersectLine.Getdirection());
    p12 = IntersectLine.GetPoint() + IntersectLine.Getdirection() * ( (tr1.getVertice(1) - IntersectLine.GetPoint()) * IntersectLine.Getdirection()) /
                                     (IntersectLine.Getdirection() * IntersectLine.Getdirection());
    p13 = IntersectLine.GetPoint() + IntersectLine.Getdirection() * ( (tr1.getVertice(2) - IntersectLine.GetPoint()) * IntersectLine.Getdirection()) /
                                     (IntersectLine.Getdirection() * IntersectLine.Getdirection());
    p21 = IntersectLine.GetPoint() + IntersectLine.Getdirection() * ( (tr2.getVertice(0) - IntersectLine.GetPoint()) * IntersectLine.Getdirection()) /
                                     (IntersectLine.Getdirection() * IntersectLine.Getdirection());
    p22 = IntersectLine.GetPoint() + IntersectLine.Getdirection() * ( (tr2.getVertice(1) - IntersectLine.GetPoint()) * IntersectLine.Getdirection()) /
                                     (IntersectLine.Getdirection() * IntersectLine.Getdirection());
    p23 = IntersectLine.GetPoint() + IntersectLine.Getdirection() * ( (tr2.getVertice(2) - IntersectLine.GetPoint()) * IntersectLine.Getdirection()) /
                                     (IntersectLine.Getdirection() * IntersectLine.Getdirection());
    if ( IntersectLine.IsOverlap(tr1.getVertice(1), tr1.getVertice(2),tr1.getVertice(3),
                                 tr2.getVertice(1),tr2.getVertice(2),tr2.getVertice(3) ))   return true;
    return false;
}

int triangle_geometry::triangle_handler::sign_of_dist(const Plane &plane, const Point &point) const {
    double dist = plane.get_norm_vec() * point;
    if ( dist < 0) return -1;
    if ( dist > 0) return 1;
    else return 0;
}

triangle_geometry::Plane::Plane() = default;

triangle_geometry::Plane::Plane(const Point& norm_vec, const Point &point){
    _norm_vec = norm_vec;
    _point = point;
    _d = norm_vec * point;
}

triangle_geometry::Plane::Plane(const triangle_geometry::triangle &tr): Plane( ((tr.getVertice(1) - tr.getVertice(2)) % (tr.getVertice(2) - tr.getVertice(3)))
                                                                               , tr.getVertice(3) ) {}

triangle_geometry::Point triangle_geometry::Plane::get_norm_vec() const{
    return _norm_vec;
}

double triangle_geometry::Plane::get_d() const{
    return _d;
}

triangle_geometry::triangle::triangle() = default;

bool triangle_geometry::triangle::is_degenerate() const{//true if degenerate
    Point norm_vec = ( _p1 - _p2, _p2 - _p3 );
    if ( norm_vec.get_x() == 0 && norm_vec.get_y() == 0 && norm_vec.get_z() == 0) return true;
    return false;
}

triangle_geometry::Point::Point(const double &x, const double &y, const double &z): _x(x), _y(y), _z(z) {}

triangle_geometry::Point::Point():Point(0.0, 0.0, 0.0)  {}

triangle_geometry::Point::Point(const Point &source) = default;

triangle_geometry::Point triangle_geometry::Point::operator-(const Point &source) const {
    return Point(_x - (source._x), _y - (source._y), _z - (source._z));
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

bool triangle_geometry::triangle_handler::is_coincident(const triangle_geometry::Plane &p1,
                                                        const triangle_geometry::Plane &p2,
                                                        const triangle_geometry::Point &point1,
                                                        const triangle_geometry::Point &point2) const {
    if ( (point1 - point2) * p1.get_norm_vec() == 0 && p1.get_d() == p2.get_d() ) return true;
    return false;
}

double triangle_geometry::Point::operator*(const Point &source) const{
    return ( get_x() * source.get_x() + get_y() * source.get_y() + get_z() * source.get_z());
}

triangle_geometry::Point triangle_geometry::Point::operator%(const Point &source) const{
    return Point(get_y() * source.get_z() - get_z() * source.get_y(), get_z() * source.get_x() - get_x() * source.get_z(),
                 get_x() * source.get_y() - get_y() * source.get_x());
}

triangle_geometry::Point triangle_geometry::Point::operator*(const int &source) const {
    return triangle_geometry::Point(_x * source, _y * source, _z * source);
}

double triangle_geometry::Point::length() const {
    return sqrt(_x * _x + _y * _y + _z * _z);
}

triangle_geometry::Point triangle_geometry::Point::operator*(const double &source) const {
    return triangle_geometry::Point(_x * source, _y * source, _z * source);
}

triangle_geometry::Point triangle_geometry::Point::operator+(const triangle_geometry::Point &source) const {
    return triangle_geometry::Point(_x + source._x, _y + source._y, _z + source._z);
}

triangle_geometry::Point triangle_geometry::Point::operator/(const double &source) const {
    return Point(_x / source, _y / source, _z / source);
}

bool triangle_geometry::triangle_handler::is_intersect2D(const triangle &tr1, const triangle &tr2) const {
    for ( int i0 = 0, i1 = 2; i0 < 3; i1 = 0, i0++){
        Point point = GetNormVector(tr1, tr1.getVertice(i0), tr1.getVertice(i1));
        int min = GetExtremeIndex(tr2, point * (-1));
        Point diff = tr2.getVertice(min) - tr1.getVertice(i0);
        if( point * diff  > 0 ) return false;
    }
    for( int i0 = 0, i1 = 2; i0 < 3; i1 = i0, i0++){
        Point point = GetNormVector(tr2, tr2.getVertice(i0), tr2.getVertice(i1));
        int min = GetExtremeIndex(tr1, point * (-1));
        Point diff = tr1.getVertice(min) - tr2.getVertice(i0);
        if ( point * diff > 0 ) return false;
    }
    return true;
}

int triangle_geometry::triangle_handler::GetMiddleIndex(const int &i0, const int &i1) const {
    if ( i0 < i1 ) return (i0 + i1)/2;
    return ( ((i0 + i1 + 3)/2 ) % 3);
}

int triangle_geometry::triangle_handler::GetExtremeIndex(const triangle &tr, const Point &point) const {
    int i0 = 0, i1 = 0;
    while ( true){
        int mid = GetMiddleIndex(i0, i1);
        int next = (mid + 1) % 3;
        Point edge = tr.getVertice(next) - tr.getVertice(mid);
        if ( point * edge > 0){
            if ( mid == i0) return i1;
            i0 = mid;
        }else{
            int prev = ( mid + 2) % 3;
            edge = tr.getVertice(mid) - tr.getVertice(prev);
            if (point * edge >= 0) return mid;
            i1 = mid;
        }
    }
}

triangle_geometry::Point triangle_geometry::triangle_handler::GetNormVector(const triangle_geometry::triangle& tr,
                                                                            const triangle_geometry::Point &p1,
                                                                            const triangle_geometry::Point &p2) const {
    return (  (p1 - p2) % Plane( tr).get_norm_vec() );
}

triangle_geometry::Line triangle_geometry::triangle_handler::GetLine(const triangle_geometry::Plane &p1,
                                                                     const triangle_geometry::Plane &p2) const {
    Point direction = p1.get_norm_vec() % p2.get_norm_vec();
    if ( direction.length() == 0) {
        std::cout<<"No line";
        exit(-1);
    }
    double s1, s2, a, b, n1n2dot, n1normsqr, n2normsqr;
    s1 = p1.get_d();
    s2 = p2.get_d();
    n1n2dot = p1.get_norm_vec() * p2.get_norm_vec();
    n1normsqr = p1.get_norm_vec() * p1.get_norm_vec();
    n2normsqr = p2.get_norm_vec() * p2.get_norm_vec();
    a = (s2 * n1n2dot - s1 * n2normsqr) / (n1n2dot * n1n2dot - n1normsqr * n2normsqr);
    b = (s1 * n1n2dot - s2 * n2normsqr) / (n1n2dot * n1n2dot - n1normsqr * n2normsqr);
    return Line{direction, p1.get_norm_vec() * a + p2.get_norm_vec() * b};
}

triangle_geometry::Point triangle_geometry::Line::GetPoint() {
    return _point;
}

triangle_geometry::Point triangle_geometry::Line::Getdirection() {
    return _direction;
}

triangle_geometry::Line::Line(const triangle_geometry::Point &direction, const triangle_geometry::Point &point):
_direction(direction), _point(point) {}

bool triangle_geometry::Line::IsOverlap(const triangle_geometry::Point &p1, const triangle_geometry::Point &p2,
                                        const triangle_geometry::Point &p3, const triangle_geometry::Point &p4,
                                        const triangle_geometry::Point &p5, const triangle_geometry::Point &p6) {
    _point = _point + _direction * 10;
    Point tr1Min, tr1Max, tr2Min, tr2Max;
    if ( (p1 - _point).length() < (p2 - _point).length()){//p1<p2
        if ( (p1 - _point).length() < (p3 - _point).length()){//p1<p3, p1<p2
            tr1Min = p1;
            if ( (p3 - _point).length() < ( (p2 - _point).length()))    tr1Max = p2;//p3<p2, p1<p3, p1<p2
            else tr1Max = p3;//p2<p3, p1<p3, p1<p2
        }
        else { //p3<p1,p1<p2
            tr1Min = p3;
            tr1Max = p2;
        }
    }else{//p2<p1
        if ( (p2 - _point).length() < (p3 - _point).length()){//p2<p1, p2<p3
            tr1Min = p2;
            if ( (p1 - _point).length() < (p3 - _point).length()){//p2<p1, p2<p3, p1<p3
                tr1Max = p3;
            }else{//p2<p1, p2<p3, p3<p1
                tr1Max = p1;
            }
        }
        else {//p2<p1, p3<p2
            tr1Min = p3;
            tr1Max =p1;
        }
    }
    if ( (p4 - _point).length() < (p5 - _point).length()){//p1<p2
        if ( (p4 - _point).length() < (p6 - _point).length()){//p1<p3, p1<p2
            tr1Min = p4;
            if ( (p6 - _point).length() < ( (p5 - _point).length()))    tr1Max = p5;//p3<p2, p1<p3, p1<p2
            else tr1Max = p6;//p2<p3, p1<p3, p1<p2
        }
        else { //p3<p1,p1<p2
            tr1Min = p6;
            tr1Max = p5;
        }
    }else{//p2<p1
        if ( (p5 - _point).length() < (p6 - _point).length()){//p2<p1, p2<p3
            tr1Min = p5;
            if ( (p4 - _point).length() < (p6 - _point).length()){//p2<p1, p2<p3, p1<p3
                tr1Max = p6;
            }else{//p2<p1, p2<p3, p3<p1
                tr1Max = p4;
            }
        }
        else {//p2<p1, p3<p2
            tr1Min = p6;
            tr1Max =p4;
        }
    }
    double tr1MinDist, tr1MaxDist, tr2MinDist, tr2MaxDist;
    tr1MinDist = (tr1Min - _point).length();
    tr1MaxDist = (tr1Max - _point).length();
    tr2MinDist = (tr2Min - _point).length();
    tr2MaxDist = (tr2Max - _point).length();
    if ( tr1MaxDist > tr2MaxDist && tr2MaxDist > tr1MinDist) return true;
    if ( tr2MaxDist > tr1MaxDist && tr1MaxDist > tr2MinDist) return true;
    if ( tr1MaxDist > tr2MaxDist && tr2MinDist > tr1MinDist) return true;
    if ( tr2MaxDist > tr1MaxDist && tr1MinDist > tr2MinDist) return true;
    if ( tr1MaxDist == tr2MaxDist && tr1MinDist == tr2MinDist) return true;
    if ( tr1MaxDist == tr2MinDist || tr2MaxDist == tr1MinDist) return true;
    return false;
}
