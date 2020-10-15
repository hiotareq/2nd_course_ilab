#include <iostream>
#include "triangle_space.hpp"
#include <fstream>

int main(){
        std::cout << "started. Wait" << std::endl;
        int NumberOfTriangles;
        float MinCoord = 0, MaxCoord = 0;
        std::vector<triangle_geometry::triangle> Triangles;
        triangle_space::Cube* FirstCube = new triangle_space::Cube;
        triangle_space::SpaceDivider* SpDiv = new triangle_space::SpaceDivider(8);
        std::cout << "Введите количество треугольников" << std::endl;
        std::cin >> NumberOfTriangles;
        std::cout << "Введите вершины треугольников" << std::endl;
        for ( int i = 0 ; i < NumberOfTriangles; i++){
            triangle_geometry::triangle NewTriangle;
            std::cin >> NewTriangle;
            Triangles.push_back(NewTriangle);
            if ( MinCoord > NewTriangle.get_min_coord()) MinCoord = NewTriangle.get_min_coord();
            if ( MaxCoord < NewTriangle.get_max_coord()) MaxCoord = NewTriangle.get_max_coord();
        }
        FirstCube->SetCube(MinCoord, MinCoord, MinCoord, MaxCoord - MinCoord);
        SpDiv->look_at_cube(*FirstCube);
        std::cout << SpDiv->get_number_of_intersections() << std::endl;
        return 0;
}
