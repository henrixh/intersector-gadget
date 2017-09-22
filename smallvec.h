#ifndef SMALLVEC_H
#define SMALLVEC_H

#ifdef __SSE4_2__
#include <immintrin.h>
#include <smmintrin.h>
#endif // __SSE4_2__

#include <iostream>

class vec3 {
private:

#ifdef __SSE4_2__
  __m128 d;
  inline vec3(const __m128 & n);
#endif // __SSE4_2__

#ifndef __SSE4_2__
  float d[4];
#endif // __SSE4_2__

public:
  inline vec3();
  inline vec3(const float v0, const float v1, const float v2);
  inline vec3(const float * vs);
  float inline length() const;
  float inline x() const;
  float inline y() const;
  float inline z() const;
  void inline normalize();
  void inline scale(const float s);
  void inline flip();

  vec3 inline friend operator+(const vec3 & lhs, const vec3 & rhs);
  vec3 inline friend operator-(const vec3 & lhs, const vec3 & rhs);
  vec3 inline & operator+=(const vec3 & rhs);
  vec3 inline & operator-=(const vec3 & rhs);

  friend vec3 inline cross(const vec3 & lhs, const vec3 & rhs);
  friend vec3 inline operator*(const vec3 & lhs, const vec3 & rhs);
  friend vec3 inline operator/(const vec3 & lhs, const vec3 & rhs);
  friend float inline dot(const vec3 & lhs, const vec3 & rhs);

  friend bool inline operator==(const vec3 & lhs, const vec3 & rhs);
  friend bool inline operator!=(const vec3 & lhs, const vec3 & rhs);
  friend vec3 inline normalize(const vec3 & v);
  friend vec3 inline scale(const vec3 & v, const float s);
  friend vec3 inline flip(const vec3 & v);
  friend vec3 inline max(const vec3 & a, const vec3 & b);
  friend vec3 inline min(const vec3 & a, const vec3 & b);
  friend float inline horizontal_max(const vec3 & a);
  friend float inline horizontal_min(const vec3 & a);
};


#ifndef __SSE4_2__
#include <math.h>

using namespace std;

inline vec3::vec3() {
  d[0] = 0;
  d[1] = 0;
  d[2] = 0;
  d[3] = 0;
}

inline vec3::vec3(float v0, float v1, float v2) {
  d[0] = v0;
  d[1] = v1;
  d[2] = v2;
  d[3] = 0.0;
}

inline vec3::vec3(const float * vs) {
  d[0] = vs[0];
  d[1] = vs[1];
  d[2] = vs[2];
  d[3] = 0.0;
}

float inline vec3::x() const {
  return d[0];
}

float inline vec3::y() const {
  return d[1];
}

float inline vec3::z() const {
  return d[2];
}

float inline vec3::length() const {
  return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

void inline vec3::normalize() {
  float scale = 1.0f/length();
  d[0] *= scale;
  d[1] *= scale;
  d[2] *= scale;
  d[3] *= scale;
}

void inline vec3::scale(const float s) {
  d[0] *= s;
  d[1] *= s;
  d[2] *= s;
  d[3] *= s;
}

void inline vec3::flip() {
  scale(-1.0);
}

vec3 inline operator+(const vec3 & lhs, const vec3 & rhs) {
  return vec3(lhs.d[0] + rhs.d[0], lhs.d[1] + rhs.d[1], lhs.d[2] + rhs.d[2]);
}

vec3 inline & vec3::operator+=(const vec3 & rhs) {
  d[0] += rhs.d[0];
  d[1] += rhs.d[1];
  d[2] += rhs.d[2];
  return *this;
}

vec3 inline & vec3::operator-=(const vec3 & rhs) {
  d[0] -= rhs.d[0];
  d[1] -= rhs.d[1];
  d[2] -= rhs.d[2];
  return *this;
}

vec3 inline operator-(const vec3 & lhs, const vec3 & rhs) {
  return vec3(lhs.d[0] - rhs.d[0], lhs.d[1] - rhs.d[1], lhs.d[2] - rhs.d[2]);
}

vec3 inline cross(const vec3 & lhs, const vec3 & rhs) {
  return vec3(lhs.d[1] * rhs.d[2] - lhs.d[2] * rhs.d[1],
              lhs.d[2] * rhs.d[0] - lhs.d[0] * rhs.d[2],
              lhs.d[0] * rhs.d[1] - lhs.d[1] * rhs.d[0]);
}

vec3 inline operator*(const vec3 & lhs, const vec3 & rhs) {
  return vec3(lhs.d[0] * rhs.d[0], lhs.d[1] * rhs.d[1], lhs.d[2] * rhs.d[2]);
}

vec3 inline operator/(const vec3 & lhs, const vec3 & rhs) {
  return vec3(lhs.d[0] / rhs.d[0], lhs.d[1] / rhs.d[1], lhs.d[2] / rhs.d[2]);
}

float inline dot(const vec3 & lhs, const vec3 & rhs) {
  float d[4];
  d[0] = lhs.d[0] * rhs.d[0];
  d[1] = lhs.d[1] * rhs.d[1];
  d[2] = lhs.d[2] * rhs.d[2];
  return d[0] + d[1] + d[2];
}

bool inline operator==(const vec3 & lhs, const vec3 & rhs) {
  return lhs.d[0] == rhs.d[0] && lhs.d[1] == rhs.d[1] && lhs.d[2] == rhs.d[2];
}

bool inline operator!=(const vec3 & lhs, const vec3 & rhs) {
  return !(lhs == rhs);
}

vec3 inline normalize(const vec3 & v) {
  vec3 res = v;
  res.normalize();
  return res;
}

vec3 inline scale(const vec3 & v, const float s) {
  vec3 res = v;
  res.scale(s);
  return res;
}

vec3 inline flip(const vec3 & v) {
  vec3 res = v;
  res.flip();
  return res;
}

vec3 inline min(const vec3 & a, const vec3 & b) {
  return vec3(min(a.x(), b.x()), min(a.y(), b.y()), min(a.z(), b.z()));
}

vec3 inline max(const vec3 & a, const vec3 & b) {
  return vec3(max(a.x(), b.x()), max(a.y(), b.y()), max(a.z(), b.z()));
}

float inline horizontal_min(const vec3 & a) {
  return std::min(a.x(), std::min(a.y(), a.z()));
}

float inline horizontal_max(const vec3 & a) {
  return std::max(a.x(), std::max(a.y(), a.z()));
}

#endif // __SSE4_2__


#ifdef __SSE4_2__

inline vec3::vec3() {
  d = _mm_set_ps(0,0,0,0);
}

inline vec3::vec3(const __m128 & n) {
  d = n;
}

inline vec3::vec3(const float v0, const float v1, const float v2) {
  d = _mm_set_ps(v0, v1, v2, 0.0f);
}

inline vec3::vec3(const float * vs) {
  d = _mm_set_ps(vs[0], vs[1], vs[2], 0.0f);
}


float inline vec3::x() const {
  return _mm_cvtss_f32(_mm_shuffle_ps(d, d, _MM_SHUFFLE(0, 0, 0, 3)));
}

float inline vec3::y() const {
  return _mm_cvtss_f32(_mm_shuffle_ps(d, d, _MM_SHUFFLE(0, 0, 0, 2)));
}

float inline vec3::z() const {
  return _mm_cvtss_f32(_mm_shuffle_ps(d, d, _MM_SHUFFLE(0, 0, 0, 1)));
}

float inline vec3::length() const {
  return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(d,d,0xE1)));
}

void inline vec3::normalize() {
  const __m128 scale = _mm_rsqrt_ps(_mm_dp_ps(d,d, 0xee));
  d = _mm_mul_ps(d, scale);
}

void inline vec3::scale(const float s) {
  d = _mm_mul_ps(d, _mm_set1_ps(s));
}

void inline vec3::flip() {
  scale(-1.0);
}

vec3 inline operator+(const vec3 & lhs, const vec3 & rhs) {
  return vec3(_mm_add_ps(lhs.d, rhs.d));
}

vec3 inline & vec3::operator+=(const vec3 & rhs) {
  this->d = this->d + rhs.d;
  return *this;
}

vec3 inline & vec3::operator-=(const vec3 & rhs) {
  this->d = this->d - rhs.d;
  return *this;
}

vec3 inline operator-(const vec3 & lhs, const vec3 & rhs) {
  return vec3(_mm_sub_ps(lhs.d, rhs.d));
}

vec3 inline cross(const vec3 & lhs, const vec3 & rhs) {
  // Borrowed from
  // http://fastcpp.blogspot.se/2011/04/vector-cross-product-using-sse-code.html

  return vec3(
    _mm_sub_ps(
               _mm_mul_ps(_mm_shuffle_ps(lhs.d, lhs.d, _MM_SHUFFLE(2, 1, 3, 0)),
                          _mm_shuffle_ps(rhs.d, rhs.d, _MM_SHUFFLE(1, 3, 2, 0))),

               _mm_mul_ps(_mm_shuffle_ps(lhs.d, lhs.d, _MM_SHUFFLE(1, 3, 2, 0)),
                          _mm_shuffle_ps(rhs.d, rhs.d, _MM_SHUFFLE(2, 1, 3, 0)))
               )
              );
}

vec3 inline operator*(const vec3 & lhs, const vec3 & rhs) {
  return vec3(_mm_mul_ps(lhs.d, rhs.d));
}

vec3 inline operator/(const vec3 & lhs, const vec3 & rhs) {
  return vec3(_mm_div_ps(lhs.d, rhs.d));
}

float inline dot(const vec3 & lhs, const vec3 & rhs) {
  return _mm_cvtss_f32(_mm_dp_ps(lhs.d, rhs.d, 0xE1)); // 0x71 is product of
                                                      // first 3, save result in
                                                      // first.
}

vec3 inline normalize(const vec3 & v) {
  const __m128 scale = _mm_rsqrt_ps(_mm_dp_ps(v.d,v.d, 0xee));
  return vec3(_mm_mul_ps(v.d, scale));
}

vec3 inline scale(const vec3 & v, const float s) {
  const __m128 scale = _mm_set1_ps(s);
  return vec3(_mm_mul_ps(v.d, scale));
}

vec3 inline flip(const vec3 & v) {
  __m128 scale = _mm_set1_ps(-1.0);
  return vec3(_mm_mul_ps(v.d, scale));
}

vec3 inline max(const vec3 & a, const vec3 & b) {
  return vec3(_mm_max_ps(a.d,b.d));
}

vec3 inline min(const vec3 & a, const vec3 & b) {
  return vec3(_mm_min_ps(a.d,b.d));
}

float inline horizontal_min(const vec3 & a) {
  return std::min(a.x(), std::min(a.y(), a.z()));
}

float inline horizontal_max(const vec3 & a) {
  return std::max(a.x(), std::max(a.y(), a.z()));
}


#endif // __SSE4_2__
#endif /* SMALLVEC_H */
