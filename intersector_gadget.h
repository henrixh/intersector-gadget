#include "smallvec.h"
#include <climits>
#include <vector>

namespace gadget {
  /////////////////////////////////////////////
  // Redefine some typenames for readability //
  /////////////////////////////////////////////

  typedef vec3 position_t;
  typedef vec3 direction_t;
  typedef vec3 normal_t;
  typedef vec3 uv_t;
  typedef size_t id_t;
  typedef float fp;

  // Enums are boring.
  typedef unsigned int axis_t;
  static const axis_t X_AXIS = 1;
  static const axis_t Y_AXIS = 2;
  static const axis_t Z_AXIS = 3;


  const static fp epsilon = 0.00001;
  const static fp infinity = SHRT_MAX;

  class ray {
  public:
    position_t origin;
    direction_t dir;
    ray() = default;
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

    vec3 center_point() const;
  };

  fp intersect_triangle(const triangle& t, const ray& r);
  intersection full_intersect_triangle(const triangle&t, const ray& r);

  class aabb {
  private:
    vec3 max_v, min_v;
    void merge_with(const aabb & other);
  public:
    aabb(const vec3 & max, const vec3 & min);
    aabb(const triangle & tri);
    aabb(const std::vector<triangle> & tris);
    inline fp intersect(const ray & r) const;
    inline fp intersect(const position_t & origin, const vec3 & inv_dir) const;
    axis_t longest_axis() const;
    fp mid_point_on_longest_axis() const;
  };

  inline fp aabb::intersect(const ray & r) const {
    vec3 inv_dir = vec3(1,1,1) / r.dir;
    return intersect(r.origin, inv_dir);
  }

  inline fp aabb::intersect(const position_t & origin, const vec3 & inv_dir) const {

    vec3 t1 = (max_v - origin) * inv_dir;
    vec3 t2 = (min_v - origin) * inv_dir;

    vec3 tmax = max(t1,t2);
    vec3 tmin = min(t1,t2);

    fp ftmax = horizontal_min(tmax);
    fp ftmin = horizontal_max(tmin);

    return ftmin + (ftmax < ftmin) * infinity + (ftmax < 0) * infinity;
  }

  struct bvh_node{
    aabb bound;
    uint32_t left, right;
    uint32_t tri_index;
  };

  class intersector {
  private:
    static const size_t threshold = 6;
    std::vector< std::vector<triangle> > tris;
    std::vector< bvh_node > nodes;

    void build_tree(uint32_t node_index, std::vector<triangle> & t);
    fp get_point(axis_t a, vec3 v);
    intersection leaf_intersect(const ray & r, triangle tri) const;

  public:
    intersector(std::vector<triangle> triangles);
    intersection intersect(const ray & r) const;
  };
}
