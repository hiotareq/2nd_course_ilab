#ifndef TRIANGLE_SPACE_HPP
#define TRIANGLE_SPACE_HPP

#include "triangle.hpp"

namespace triangle_space{
    class Cube{
        private:
        double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax, length;
        std::vector<triangle_geometry::triangle> triangles;
        public:
        Cube();
        Cube(const double& xmin, const double& ymin, const double& zmin, const double& xmax, const double& ymax, const double& zmax);

        std::vector<triangle_geometry::triangle> GetTriangles() const;

        void TriangleToCube(const triangle_geometry::triangle& tr, std::vector<Cube&>& cubes) const;
    };

    class SpaceDivider{
        private:
            int CounterOfIntersections;
            int threshold;
        public:
            SpaceDivider();
            void LookAtCube(const Cube& cube);
            Cube SetCube();

            void DivideAndLook( const Cube& cube);

            
    };

}


#endif