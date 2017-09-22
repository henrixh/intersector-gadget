////////////////////////////////////////////////////////////////////////////
// DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE                            //
// Version 2, December 2004                                               //
//                                                                        //
//   Copyright (C) 2004 Sam Hocevar <sam@hocevar.net>                     //
//                                                                        //
//   Everyone is permitted to copy and distribute verbatim or modified    //
//   copies of this license document, and changing it is allowed as long  //
//   as the name is changed.                                              //
//                                                                        //
//   DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE                          //
//   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION      //
//                                                                        //
//   0. You just DO WHAT THE FUCK YOU WANT TO.                            //
////////////////////////////////////////////////////////////////////////////

// https://github.com/henrixh/intersector-gadget
// Maintainer: Henrik Henriksson <hx@hx.ax>
// Authors: Mattias Hammar <mahavilan@gmail.com> & Henrik Henriksson <hx@hx.ax>

#include "intersector_gadget.h"
#include <cmath>
#include <stack>
#include <algorithm>
#include <utility>

namespace gadget {

  ray::ray(position_t o, direction_t d) {
    origin = o;
    dir = d;
  }

  fp intersect_triangle(const triangle& tri, const ray& r) {
    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    vec3 e1 = tri.v1 - tri.v0;
    vec3 e2 = tri.v2 - tri.v0;
    vec3 p = cross(r.dir, e2);
    fp det = dot(e1, p);

    if (det == 0) return infinity;

    fp inv_det = 1.0f / det;
    vec3 t = r.origin - tri.v0;
    fp u = dot(t, p) * inv_det;

    if (u < 0.0f || u > 1.0f) return infinity;


    vec3 q = cross(t, e1);
    fp v = dot(r.dir, q) * inv_det;

    if (v < 0.0f || u + v > 1.0f) return infinity;


    fp dist = dot(e2, q) * inv_det;
    if (dist > epsilon) {
      return dist - epsilon;
    }

    return infinity;
  }


  intersection full_intersect_triangle(const triangle&tri, const ray& r) {
    fp dist = intersect_triangle(tri, r);

    intersection res;
    if (dist == infinity) {
      res.hit = false;
      return res;
    }
    res.hit = true;
    res.hit_position = scale(r.dir, dist) + r.origin;

    // http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    vec3 v0 = tri.v1 - tri.v0;
    vec3 v1 = tri.v2 - tri.v0;
    vec3 v2 = res.hit_position - tri.v0;

    fp d00 = dot(v0,v0);
    fp d01 = dot(v0,v1);
    fp d11 = dot(v1,v1);
    fp d20 = dot(v2,v0);
    fp d21 = dot(v2,v1);
    fp denom = d00 * d11 - d01 * d01;
    fp v = (d11 * d20 - d01 * d21) / denom;
    fp w = (d00 * d21 - d01 * d20) / denom;
    fp u = 1.0 - v - w;

    res.hit_normal = normalize(scale(tri.n0,u) + scale(tri.n1, v) + scale(tri.n2, w));
    res.hit_normal = normalize(cross(tri.v1 - tri.v0, tri.v2 - tri.v0));
    res.uv = scale(tri.t0,u) + scale(tri.t1, v) + scale(tri.t2, w);
    res.incoming_ray = r;

    return res;
  }

  vec3 triangle::center_point() const {
    return scale(v0 + v1 + v2, 1.0/3.0);
  }

aabb::aabb(const vec3 & max_v, const vec3 & min_v) {
  this->max_v = max_v;
  this->min_v = min_v;
}

aabb::aabb(const triangle & tri) {
  fp max_x = std::max(tri.v0.x(), std::max(tri.v1.x(), tri.v2.x()));
  fp max_y = std::max(tri.v0.y(), std::max(tri.v1.y(), tri.v2.y()));
  fp max_z = std::max(tri.v0.z(), std::max(tri.v1.z(), tri.v2.z()));
  fp min_x = std::min(tri.v0.x(), std::min(tri.v1.x(), tri.v2.x()));
  fp min_y = std::min(tri.v0.y(), std::min(tri.v1.y(), tri.v2.y()));
  fp min_z = std::min(tri.v0.z(), std::min(tri.v1.z(), tri.v2.z()));

  max_v = {max_x, max_y, max_z};
  min_v = {min_x, min_y, min_z};
}

  aabb::aabb(const std::vector<triangle> & tris) {
  if (!tris.size()) {
    return;
  }

  aabb res {tris[0]};
  for (size_t i = 1; i < tris.size(); ++i) {
    res.merge_with({tris[i]});
  }

  max_v = res.max_v;
  min_v = res.min_v;
}

void aabb::merge_with(const aabb & other) {
  fp max_x = std::max(other.max_v.x(), max_v.x());
  fp max_y = std::max(other.max_v.y(), max_v.y());
  fp max_z = std::max(other.max_v.z(), max_v.z());
  fp min_x = std::min(other.min_v.x(), min_v.x());
  fp min_y = std::min(other.min_v.y(), min_v.y());
  fp min_z = std::min(other.min_v.z(), min_v.z());

  max_v = {max_x, max_y, max_z};
  min_v = {min_x, min_y, min_z};
}


axis_t aabb::longest_axis() const {
  vec3 diff = max_v - min_v;
  if (diff.x() >= diff.y() && diff.x() >= diff.z()) return X_AXIS;
  if (diff.y() >= diff.z() && diff.y() >= diff.x()) return Y_AXIS;
  if (diff.z() >= diff.x() && diff.z() >= diff.y()) return Z_AXIS;
  return X_AXIS;
}

fp aabb::mid_point_on_longest_axis() const {
  axis_t axis = longest_axis();

  if (axis == X_AXIS) {
    return min_v.x() + 0.5*(max_v.x() - min_v.x());
  }
  if (axis == Y_AXIS) {
    return min_v.y() + 0.5*(max_v.y() - min_v.y());
  }
  if (axis == Z_AXIS) {
    return min_v.z() + 0.5*(max_v.z() - min_v.z());
  }

  // Should never happen.
  return infinity;
}

  intersector::intersector(std::vector<triangle> triangles){
    std::cout << "Creating bvh tree with " << triangles.size() << " triangles" << std::endl;
  bvh_node root {triangles, 0, 0, 0};
  nodes.push_back(root);
  build_tree(0, triangles);
  std::cout << "Created bvh tree with " << nodes.size() << " nodes" << std::endl;
}

fp intersector::get_point(axis_t a, vec3 v){
  if (a == X_AXIS){
    return v.x();
  }else if (a == Y_AXIS){
    return v.y();
  }else{
    return v.z();
  }
}

  void intersector::build_tree(uint32_t node_index, std::vector<triangle> & t){
  bvh_node node = nodes[node_index];
  if(t.size() <= threshold){
    node.tri_index = tris.size();
    tris.push_back(t);
    nodes[node_index] = node;
    return;
  }

  axis_t axis = node.bound.longest_axis();
  float mid_point = node.bound.mid_point_on_longest_axis();


  std::vector<triangle> left_vec;
  std::vector<triangle> right_vec;

  for (auto && a : t) {
    if (get_point(axis, a.center_point()) < mid_point) {
      left_vec.push_back(a);
    } else {
      right_vec.push_back(a);
    }
  }

  if (right_vec.empty() || left_vec.empty()) {
    left_vec.clear(); left_vec.shrink_to_fit();
    right_vec.clear(); right_vec.shrink_to_fit();
    std::sort(t.begin(), t.end(), [&](const triangle & a, const triangle & b) -> bool {
        return get_point(axis, a.center_point()) < get_point(axis, b.center_point());
      });
    left_vec = {t.begin(), t.begin() + t.size()/2};
    right_vec = {t.begin() + t.size()/2, t.end()};
  }

  t.clear(); t.shrink_to_fit();

  node.left = nodes.size();
  bvh_node left_node  {left_vec,  0, 0, 0};
  nodes.push_back(left_node);

  node.right = nodes.size();
  bvh_node right_node {right_vec, 0, 0, 0};
  nodes.push_back(right_node);

  build_tree(node.left, left_vec);
  build_tree(node.right, right_vec);
  nodes[node_index] = node;
}

intersection intersector::leaf_intersect(const ray & r, triangle tri) const {
  intersection res;
  res = full_intersect_triangle(tri, r);
  res.incoming_ray = r;
  res.id = tri.id;
  res.hit = true;

  return res;
}

intersection intersector::intersect(const ray & r) const {
  std::stack<uint32_t> frontier;
  vec3 inv_ray_dir = vec3(1,1,1) / r.dir;
  frontier.push(0);

  fp closest_distance = infinity;
  std::pair<size_t, size_t> closest_object {};


  while (frontier.size()) {
    uint32_t node_index = frontier.top(); frontier.pop();
    bvh_node n = nodes[node_index];

    // 0 in any pointer marks a leaf
    if (n.left == 0 || n.right == 0){
      if (n.bound.intersect(r.origin, inv_ray_dir) < closest_distance) {
        for (size_t i { 0 }; i < tris[n.tri_index].size(); ++i) {
          fp hit = intersect_triangle(tris[n.tri_index][i], r);
          if (hit < closest_distance) {
            closest_distance = hit;
            closest_object = {n.tri_index, i};
          }
        }
      }
      continue;
    }

    fp left_dist = nodes[n.left].bound.intersect(r.origin, inv_ray_dir);
    fp right_dist = nodes[n.right].bound.intersect(r.origin, inv_ray_dir);

    if (left_dist < right_dist) {
      if (right_dist < infinity) { frontier.push(n.right); }
      if (left_dist < infinity) { frontier.push(n.left); }
    } else {
      if (left_dist < infinity) { frontier.push(n.left); }
      if (right_dist < infinity) { frontier.push(n.right); }
    }
  }

  intersection tri_intersection {};
  if (closest_distance < infinity){
    tri_intersection = leaf_intersect(r, tris[closest_object.first][closest_object.second]);
  }
  return tri_intersection;
}

}
