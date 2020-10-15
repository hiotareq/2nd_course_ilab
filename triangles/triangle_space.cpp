#include "trianglespace.hpp"

std::vector<triangle_geometry::triangle> triangle_space::Cube::GetTriangles() const{
    return triangles;
}

void triangle_space::Cube::TriangleToCube(const triangle_geometry::triangle& tr, std::vector<Cube>& cubes) const{
    int imax = 1, jmax = 1, kmax = 1;

    if ((tr.get_max_x() - tr.get_min_x())/(_length/2) == 1) imax = 2;
    if ((tr.get_max_y() - tr.get_min_y())/(_length/2) == 1) jmax = 2;
    if ((tr.get_max_z() - tr.get_min_z())/(_length/2) == 1) kmax = 2;

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){
            for (int k = 0; k < kmax; k++){
                cubes[i + j +k].triangles.push_back(tr);
            }   
        }
    }
}

void triangle_space::Cube::SetCube(const float& xmin, const float& ymin, const float& zmin, const float& length){
    _xmin = xmin;
    _ymin = ymin;
    _zmin = zmin;
    _length = length;
}

triangle_space::SpaceDivider::SpaceDivider(const int& n) : threshold(n), CounterOfIntersections(0){}

void triangle_space::SpaceDivider::divide_and_look(const Cube& cube){
    std::vector<triangle_space::Cube> cubes(8);

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; i++){
            for (int k = 0; k < 2; k++){
                cubes[i + j + k].SetCube(cube.get_min_x() + i * cube.get_length() / 2,cube.get_min_y() + j * cube.get_length()/2, cube.get_min_z() + k * cube.get_length()/2, cube.get_length()/2);
            }
        }
    }

    for ( auto it = cube.GetTriangles().begin(); it != cube.GetTriangles().end(); it++){
        cube.TriangleToCube(*it, cubes);
    }
    
    for (int i = 0; i < 2; i++){
        for (int j = 0; j< 2; j++){
            for (int k = 0; k < 2; k++){
                look_at_cube(cubes[i + j + k]);
            }  
        }
    }
}

void triangle_space::SpaceDivider::look_at_cube(const Cube& cube){
    if (cube.GetTriangles().size() <= threshold){
        for ( auto it1 = cube.GetTriangles().begin(); it1 != cube.GetTriangles().end(); it1++){
            for ( auto it2 = cube.GetTriangles().begin(); it2 != cube.GetTriangles().end(); it2++){
                if (it1 == it2) continue;
                else{
                    bool check = false;
                    check = triangle_geometry::is_intersect(*it1, *it2);
                    if (check) CounterOfIntersections++;
                    else continue;
                }
            }
        }
    }else{
        divide_and_look(cube);
    }
}

triangle_space::Cube::Cube(const float& xmin, const float& ymin,const float& zmin, const float& length):
_xmin(xmin), _ymin(ymin), _zmin(zmin), _length(length){}

triangle_space::Cube::Cube() : Cube(0.0 ,0.0 ,0.0 ,0.0){}

int triangle_space::SpaceDivider::get_number_of_intersections() const{
    return CounterOfIntersections;
}

float triangle_space::Cube::get_min_x() const{
    return _xmin;
}

float triangle_space::Cube::get_min_y() const{
    return _ymin;
}

float triangle_space::Cube::get_min_z() const{
    return _zmin;
}

float triangle_space::Cube::get_length() const{
    return _length;
}