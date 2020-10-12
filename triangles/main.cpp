#include "triangle.hpp"
#include <iostream>

int main(){
    triangle_geometry::Point p00(0, 0, 0), p01(1, 0, 0), p02(0, 0, 1), p10(1, 1, 1), p11(0, 1,0 ), p12(0, 0, 1);
    triangle_geometry::triangle tr0(p00, p01, p02), tr1(p10, p11, p12);
    if ( triangle_geometry::is_intersect(tr0, tr1)) std::cout << "true";
    else std::cout<<"false";
    return 0;
}