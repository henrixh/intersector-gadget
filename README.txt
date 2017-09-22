This is a small, simple and somewhat decent ray/mesh intersection library. If
something more production worthy is needed, look somewhere else.

The goal is not to be correct, not to be fast. 

The goal is to be easy to use when writing toy raytracers, pathtracers or other
things where some ray/mesh intersection is needed.



Usage:

Move all .cpp and .h files to your project and add to your build system. If you
are building on x86-64, you should add '-march=native' to your flags and use the
intrinsic version of 'smallvec'. If this isn't working, remove the line
> #define SMALLVEC_USE_INTRINSIC'
from 'smallvec.h'.

Compiles with '$ g++ -O3 -pedantic -Wall -Wextra -Werror -march=native -std=gnu++14'




API:

namespace gadget {
  typedef vec3 position_t;
  typedef vec3 direction_t;
  typedef vec3 normal_t;
  typedef vec3 uv_t;
  typedef size_t id_t;
  typedef float fp;

  class ray {
  public:
    position_t origin;
    direction_t dir;
    ray(position_t o, direction_t d);
  };

  struct intersection {
    normal_t hit_normal;
    position_t hit_position;
    ray incoming_ray;
    uv_t uv;
    id_t id;
    bool hit;
  };

  class triangle {
  public:
    position_t v0, v1, v2;
    normal_t n0, n1, n2;
    uv_t t0, t1, t2;
    id_t id;
  };

  class intersector {
  public:
    intersector(std::vector<triangle> triangles);
    intersection intersect(const ray & r) const;
  };
}




Smallvec:

'smallvec.h' is also available at https://github.com/henrixh/smallvec as a
separate project, when only some 3D vectors are needed



Written by:

Mattias Hammar <mahavilan@gmail.com>
Henrik Henriksson <hx@hx.ax>
