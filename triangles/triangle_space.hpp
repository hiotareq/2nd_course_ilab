#ifndef TRIANGLE_SPACE_HPP
#define TRIANGLE_SPACE_HPP

#include "triangle.hpp"

namespace triangle_space{
    class Cube{
        private:
        float _xmin, _ymin, _zmin, _length;
        std::vector<triangle_geometry::triangle> triangles;
        public:
        Cube();
        Cube(const float& xmin, const float& ymin, const float& zmin, const float& length);
        void SetCube(const float& xmin, const float& ymin, const float& zmin, const float& length);

        float get_min_x() const;
        float get_min_y() const;
        float get_min_z() const;

        float get_length() const;

        std::vector<triangle_geometry::triangle> GetTriangles() const;

        void TriangleToCube(const triangle_geometry::triangle& tr, std::vector<Cube&>& cubes) const;
    };

    class SpaceDivider{
        private:
            int CounterOfIntersections;
            int threshold;
        public:
            SpaceDivider();
            void look_at_cube(const Cube& cube);
            
            void divide_and_look( const Cube& cube);            
    };
}


#endif