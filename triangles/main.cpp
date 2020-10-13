#include "triangle.hpp"
#include <iostream>

int main(){
    triangle_geometry::Point p00(0, 0, 0), p01(10, 0, 0), p02(0, 10, 0), p10(2, 2, 0), p11(3, 4, 0), p12(4, 2, 0);
    triangle_geometry::triangle tr0(p00, p01, p02), tr1(p10, p11, p12);
    std::boolalpha;
    std::cout << triangle_geometry::is_intersect(tr0, tr1);
    return 0;
}