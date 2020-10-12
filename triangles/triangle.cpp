#include "triangle.hpp"
#include <cmath>
#include <algorithm>

bool triangle_geometry::is_intersect(const triangle &tr1, const triangle &tr2) {
    if ( tr1.is_degenerate() ){
        /*
         * удаление вырожденного треугольника из базы данных
         */
        return false;
    }
    if ( tr2.is_degenerate()){
        /*
         * удаление вырожденного треугольника из базы данных
         */
        return false;
    }

    Plane plane_1(tr1);

    int sign11, sign12, sign13;
    sign11 = sign_of_dist(plane_1, tr2.getVertice(1));
    sign12 = sign_of_dist(plane_1, tr2.getVertice(2));
    sign13 = sign_of_dist(plane_1, tr2.getVertice(3));
    if ( sign11 * sign12 > 0 && sign11 * sign13 > 0) return false;//треугольник полностью с одной стороны от плоскости
    //else continue
    
    Plane plane_2 = Plane( tr2);
    if ( is_coincident( plane_1, plane_2, tr1.getVertice(1), tr2.getVertice(1) )) {
        if (plane_1.get_d() == plane_2.get_d()) {
            if ( is_intersect2D(tr1, tr2)) return true;
            return false;
        }return false;//треугольники в параллельных плоскостях
    }//else continue

    int sign21, sign22, sign23;
    sign21 = sign_of_dist(plane_2, tr1.getVertice(1));
    sign22 = sign_of_dist(plane_2, tr1.getVertice(2));
    sign23 = sign_of_dist(plane_2, tr1.getVertice(3));
    if ( sign21 * sign22 > 0 && sign21 * sign23 > 0) return false;//треугольник полностью с одной стороны от плоскости
    //else continue

    Vector3D edge00, edge01, edge10, edge11;
    Point p00, p01, p10, p11;//это точки для вычисления отрезков-проекций на линию пересечения

    p00 = tr1.getVertice(1);
    p01 = tr1.getVertice(3);
    if ( sign11 * sign12 < 0){//1 и 2 по разные стороны от плоскости
        edge00 = Vector3D(tr1.getVertice(1), tr1.getVertice(2));
        if ( sign11 * sign13 < 0){//1 и 3 по разные стороны от плоскости
            edge01 = Vector3D(tr1.getVertice(1), tr1.getVertice(3));
        }else{//1 и 3 по одну сторону от плоскости, с другой стороны - 2
            edge01 = Vector3D(tr1.getVertice(2), tr1.getVertice(3));
        }
    }else{//1 и 2 по одну сторону от плоскости
        edge00 = Vector3D(tr1.getVertice(1), tr1.getVertice(3));
        edge01 = Vector3D(tr1.getVertice(2), tr1.getVertice(3));
    }

    p10 = tr2.getVertice(1);
    p11 = tr2.getVertice(3);
    if ( sign21 * sign22 < 0){//1 и 2 по разные стороны от плоскости
        edge00 = Vector3D(tr2.getVertice(1), tr2.getVertice(2));
        if ( sign21 * sign23 < 0){//1 и 3 по разные стороны от плоскости
            edge11 = Vector3D(tr2.getVertice(1), tr2.getVertice(3));
        }else{//1 и 3 по одну сторону от плоскости, с другой стороны - 2
            edge11 = Vector3D(tr2.getVertice(2), tr2.getVertice(3));
        }
    }else{//1 и 2 по одну сторону от плоскости
        edge10 = Vector3D(tr1.getVertice(1), tr1.getVertice(3));
        edge11 = Vector3D(tr1.getVertice(2), tr1.getVertice(3));
    }

    //проекция точек на линию пересечения
    Line IntersectLine = GetLine(plane_1, plane_2);
    
    Point t00, t01, t10, t11;//t0 - первый треугольник и его точки пересечения с прямой, t1 - второй треугольник
    t00 = IntersectionEdgeLine(p00, IntersectLine, edge00);
    t01 = IntersectionEdgeLine(p01, IntersectLine, edge01);
    t10 = IntersectionEdgeLine(p10, IntersectLine, edge10);
    t11 = IntersectionEdgeLine(p11, IntersectLine, edge11);
    std::vector<triangle_geometry::Point> points = triangle_geometry::DefinePoints(t00, t01, t10, t11);//вектор с упорядоченными вершинами
    
    if ( (points[0] == t00 || points[0] == t01) && (points[3] == t00 || points[3] == t01)) return true;
    if ( (points[0] == t10 || points[0] == t11) && (points[3] == t10 || points[3] == t11)) return true;
    if ( (points[0] == t00 || points[0] == t01) && (points[1] == t10 || points[1] == t11)) return true;
    if ( (points[0] == t10 || points[0] == t11) && (points[1] == t00 || points[1] == t01)) return true;
    if ( t00 == t01 || t00 == t11 || t10 == t01 || t01 == t11) return true;
    return false;
}

int triangle_geometry::sign_of_dist(const Plane &plane, const Point &point) {
    float dist = plane.get_normal() * point;
    if (dist < 0) return -1;
    if (dist > 0) return 1;
    else return 0;
}

triangle_geometry::Plane::Plane() = default;

triangle_geometry::Plane::Plane(const Vector3D &norm_vec, const Point &point) {
    normal = norm_vec;
    _point = point;
    _d = (norm_vec * point)/(1e-5);
}

triangle_geometry::Plane::Plane(const triangle_geometry::triangle &tr) {
            Vector3D v1(tr.getVertice(1), tr.getVertice(2)), v2(tr.getVertice(2), tr.getVertice(3));
            normal = v1 % v2;
            _point = tr.getVertice(3);
            _d = v1 * normal;
        }

triangle_geometry::Vector3D triangle_geometry::Plane::get_normal() const {
    return normal;
}

float triangle_geometry::Plane::get_d() const {
    return _d;
}

triangle_geometry::triangle::triangle() = default;

bool triangle_geometry::triangle::is_degenerate() const {//true if degenerate
    Vector3D norm_vec(0,0,0), v1(_p1, _p2), v2(_p2, _p3);
    norm_vec = v1 % v2;
    if (norm_vec.getX() == 0 && norm_vec.getY() == 0 && norm_vec.getZ() == 0) return true;
    return false;
}

triangle_geometry::Point::Point(const float &x, const float &y, const float &z) : _x(x), _y(y), _z(z) {}

triangle_geometry::Point::Point() : Point(0.0, 0.0, 0.0) {}

float triangle_geometry::Point::get_x() const {
    return _x;
}

float triangle_geometry::Point::get_y() const {
    return _y;
}

float triangle_geometry::Point::get_z() const {
    return _z;
}

triangle_geometry::Point triangle_geometry::triangle::getVertice(const int &i) const {
    if (i == 1) return _p1;
    if (i == 2) return _p2;
    return _p3;
}

triangle_geometry::triangle::triangle(const triangle_geometry::Point &p1, const triangle_geometry::Point &p2,
                                      const triangle_geometry::Point &p3) : _p1(p1), _p2(p2), _p3(p3) {}

bool triangle_geometry::is_coincident(const triangle_geometry::Plane &p1,
                                                        const triangle_geometry::Plane &p2,
                                                        const triangle_geometry::Point &point1,
                                                        const triangle_geometry::Point &point2) {
                                                            Vector3D v(point1, point2);
    if (v * p1.get_normal() == 0 && p1.get_d() == p2.get_d()) return true;
    return false;
}


bool triangle_geometry::is_intersect2D(const triangle &tr1, const triangle &tr2) {
    for (int i0 = 0, i1 = 2; i0 < 3; i1 = 0, i0++) {
        Vector3D point = GetNormVector(tr1, tr1.getVertice(i0), tr1.getVertice(i1));
        int min = GetExtremeIndex(tr2, point * (-1));
        Vector3D diff(tr2.getVertice(min), tr1.getVertice(i0));
        if (point * diff > 0) return false;
    }
    for (int i0 = 0, i1 = 2; i0 < 3; i1 = i0, i0++) {
        Vector3D point = GetNormVector(tr2, tr2.getVertice(i0), tr2.getVertice(i1));
        int min = GetExtremeIndex(tr1, point * (-1));
        Vector3D diff(tr1.getVertice(min), tr2.getVertice(i0));
        if (point * diff > 0) return false;
    }
    return true;
}

int triangle_geometry::GetMiddleIndex(const int &i0, const int &i1) {
    if (i0 < i1) return (i0 + i1) / 2;
    return (((i0 + i1 + 3) / 2) % 3);
}

int triangle_geometry::GetExtremeIndex(const triangle &tr, const Vector3D &point) {
    int i0 = 0, i1 = 0;
    while (true) {
        int mid = GetMiddleIndex(i0, i1);
        int next = (mid + 1) % 3;
        Vector3D edge{tr.getVertice(next), tr.getVertice(mid)};
        if (point * edge > 0) {
            if (mid == i0) return i1;
            i0 = mid;
        } else {
            int prev = (mid + 2) % 3;
            edge = Vector3D(tr.getVertice(mid), tr.getVertice(prev));
            if (point * edge >= 0) return mid;
            i1 = mid;
        }
    }
}

triangle_geometry::Vector3D triangle_geometry::GetNormVector(const triangle_geometry::triangle &tr,
const triangle_geometry::Point &p1,const triangle_geometry::Point &p2) {
    Vector3D v1(p1), v2(p2);
    return ((v1 - v2) % Plane(tr).get_normal());
}

triangle_geometry::Line triangle_geometry::GetLine(const triangle_geometry::Plane &p1,
                                                                     const triangle_geometry::Plane &p2) {                                                                         
    Vector3D direction = p1.get_normal() % p2.get_normal();
    if (direction.length() == 0) {
        std::cout << "No line";
        exit(-1);
    }
    float s1, s2, a, b, n1n2dot, n1normsqr, n2normsqr;
    s1 = p1.get_d();
    s2 = p2.get_d();
    n1n2dot = p1.get_normal() * p2.get_normal();
    n1normsqr = p1.get_normal() * p1.get_normal();
    n2normsqr = p2.get_normal() * p2.get_normal();
    a = (s2 * n1n2dot - s1 * n2normsqr) / (n1n2dot * n1n2dot - n1normsqr * n2normsqr);
    b = (s1 * n1n2dot - s2 * n2normsqr) / (n1n2dot * n1n2dot - n1normsqr * n2normsqr);
    Point p = MakePointFromVector(p1.get_normal() * a + p2.get_normal() * b);
    return Line(direction, p) ;
}

triangle_geometry::Point triangle_geometry::Line::GetPoint() const{
    return _point;
}

triangle_geometry::Vector3D triangle_geometry::Line::Getdirection() const{
    return _direction;
}

triangle_geometry::Line::Line(const triangle_geometry::Vector3D &direction, const triangle_geometry::Point &point) :
        _direction(direction), _point(point) {}

triangle_geometry::Point triangle_geometry::IntersectionEdgeLine(const Point& PointFromTriangle,const Line& line,
        const Vector3D& EdgeDir) {
    Vector3D v(0, -EdgeDir.getZ(), EdgeDir.getY()), w( line.GetPoint(), PointFromTriangle);
    float s = ((v % w).length())/((v % line.Getdirection()).length());
    return ( line.GetPoint() + s * line.Getdirection());
}

triangle_geometry::Vector3D::Vector3D() : Vector3D(0, 0, 0) {
}

triangle_geometry::Vector3D::Vector3D(float _x, float _y, float _z) : _x(_x), _y(_y), _z(_z) {
}

triangle_geometry::Vector3D::~Vector3D() {
}

float triangle_geometry::Vector3D::getX() const {
    return _x;
}

float triangle_geometry::Vector3D::getY() const {
    return _y;
}

float triangle_geometry::Vector3D::getZ() const {
    return _z;
}

float triangle_geometry::Vector3D::length() const {
    return sqrt(_x * _x + _y * _y + _z * _z);
}

void triangle_geometry::Vector3D::setX(float x) {
    _x = x;
}

void triangle_geometry::Vector3D::setY(float y) {
    _y = y;
}

void triangle_geometry::Vector3D::setZ(float z) {
    _z = z;
}

bool triangle_geometry::Vector3D::operator==(const triangle_geometry::Vector3D &v) const {
    float E = 0.000001;
    return (abs(v.getX() - _x) < E) && (abs(v.getY() - _y) < E) && (abs(v.getZ() - _z) < E);
}

bool triangle_geometry::Vector3D::operator!=(const triangle_geometry::Vector3D &v) const {
    return !(*this == v);
}

triangle_geometry::Vector3D triangle_geometry::Vector3D::operator+(const triangle_geometry::Vector3D &v) const {
    return Vector3D(_x + v.getX(), _y + v.getY(), _z + v.getZ());
}

triangle_geometry::Vector3D triangle_geometry::Vector3D::operator-(const triangle_geometry::Vector3D &v) const {
    return Vector3D(_x - v.getX(), _y - v.getY(), _z - v.getZ());
}

triangle_geometry::Vector3D triangle_geometry::Vector3D::operator*(const int a) const {
    return Vector3D(_x * a, _y * a, _z * a);
}

float triangle_geometry::Vector3D::operator*(const triangle_geometry::Vector3D &v) const {
    return (_x * v.getX() + _y * v.getY() + _z * v.getZ());
}


triangle_geometry::Vector3D operator*(const int& a, const triangle_geometry::Vector3D& v) {
    return triangle_geometry::Vector3D(v.getX() * a, v.getY() * a, v.getZ() * a);
}

triangle_geometry::Vector3D operator*(const float& a, const triangle_geometry::Vector3D& v){
    return triangle_geometry::Vector3D(v.getX() * a, v.getY() * a, v.getZ() * a);
}

triangle_geometry::Vector3D triangle_geometry::Vector3D::operator*(const float a) const {
    return Vector3D(_x * a, _y * a, _z * a);
}

std::ostream &operator<<(std::ostream &os, const triangle_geometry::Vector3D &v) {
    os << "(" << v.getX() << "; " << v.getY() << "; " << v.getZ() << ")";
    return os;
}

std::istream &operator>>(std::istream &is, triangle_geometry::Vector3D &v) {
    float x, y, z;
    is >> x >> y >> z;
    v.setX(x);
    v.setY(y);
    v.setZ(z);
    return is;
}

triangle_geometry::Vector3D::Vector3D(const triangle_geometry::Point& point) : _x(point.get_x()), _y(point.get_y()), _z(point.get_z()) {}

triangle_geometry::Vector3D::Vector3D(const triangle_geometry::Point &point1, const triangle_geometry::Point &point2) :
        _x(point1.get_x() - point2.get_x()), _y(point1.get_y() - point2.get_y()), _z(point1.get_z() - point2.get_z()) {}

triangle_geometry::Vector3D triangle_geometry::Vector3D::operator%(const Vector3D& v) const {
    return Vector3D( (_y * v.getZ() - _z * v.getY())*(1000)/(1000), (_z * v.getX() - _x * v.getZ())*(1000)/(1000), (_x * v.getY() - _y * v.getX())*(1000)/(1000));
}

triangle_geometry::Point triangle_geometry::operator+(const Point& p, const Vector3D& v ){
    return Point( v.getX() + p.get_x(), v.getY() + p.get_y(), v.getZ() + p.get_z());
}

triangle_geometry::Point triangle_geometry::MakePointFromVector(const Vector3D& v){
    return Point(v.getX(), v.getY(), v.getZ());
}

triangle_geometry::Vector3D triangle_geometry::MakeVectorFromPoint(const Point& p){
    return Vector3D(p.get_x(), p.get_y(), p.get_z());
}

bool triangle_geometry::Point::operator==(const Point& p) const{
    return (_x == p.get_x() && _y == p.get_y() && _z == p.get_z());
}

std::vector<triangle_geometry::Point> triangle_geometry::DefinePoints(const Point& p1, const Point& p2, const Point& p3, const Point& p4){
    Point most_far = triangle_geometry::MostFarPoint(p1, p2, p3, p4);
    if ( most_far == p2){//2 дальше всех от 1
        return triangle_geometry::MakeVector(p2, p1, p3, p4);
    }
    if ( most_far == p3){
        return triangle_geometry::MakeVector(p3, p1, p2, p4);
    }else{//most_far == p4
        return triangle_geometry::MakeVector(p4, p1, p2, p3);
    }
}


std::vector<triangle_geometry::Point> triangle_geometry::MakeVector(const triangle_geometry::Point& most_far, 
const triangle_geometry::Point& p1, const triangle_geometry::Point& p2, const triangle_geometry::Point& p3){
    triangle_geometry::Vector3D v1(most_far, p1);
    triangle_geometry::Vector3D v2(most_far, p2), v3(most_far, p3);
    if( v2.length() < v1.length() && v3.length() < v1.length()){//most_far и 1 - крайние точки
            if ( v2.length() < v3.length()) return std::vector<triangle_geometry::Point>{p1, p3, p2, most_far};
            else return std::vector<triangle_geometry::Point>{p1, p2, p3, most_far};
    }else{
        if ( v1.length() < v2.length() && v3.length() < v2.length()){//most_far и 2 - крайние точки
            if ( v1.length() < v3.length() ) return std::vector<triangle_geometry::Point>{p2, p3, p1, most_far};
            else return std::vector<triangle_geometry::Point>{p2, p1, p3, most_far};
        }else{//most_far и 3 - крайние точки
            if ( v1.length() < v2.length() ) return std::vector<triangle_geometry::Point>{p3, p2, p1, most_far};
            else return std::vector<triangle_geometry::Point>{p3, p1, p2, most_far};
        }
    }
}

triangle_geometry::Point triangle_geometry::MostFarPoint(const triangle_geometry::Point& point_to_cmp, const triangle_geometry::Point& p1,
 const triangle_geometry::Point& p2, const triangle_geometry::Point& p3){
    triangle_geometry::Vector3D v1(point_to_cmp, p1), v2(point_to_cmp, p2), v3(point_to_cmp, p3);
    float l1 = v1.length(), l2 = v2.length(), l3 = v3.length();
    if ( l1 > l2 && l1 > l3) return p1;
    if ( l2 > l1 && l2 > l3) return p2;
    return p3;
}