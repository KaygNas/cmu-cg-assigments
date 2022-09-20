
#pragma once

#include <algorithm>
#include <cmath>
#include <ostream>

#include "log.h"

struct Vec2 {

	Vec2() {
		x = 0.0f;
		y = 0.0f;
	}
	explicit Vec2(float _x, float _y) {
		x = _x;
		y = _y;
	}
	explicit Vec2(float f) {
		x = y = f;
	}
	explicit Vec2(int32_t _x, int32_t _y) {
		x = static_cast<float>(_x);
		y = static_cast<float>(_y);
	}

	Vec2(const Vec2&) = default;
	Vec2& operator=(const Vec2&) = default;
	~Vec2() = default;

	float& operator[](uint32_t idx) {
		assert(idx <= 1);
		return data[idx];
	}
	float operator[](uint32_t idx) const {
		assert(idx <= 1);
		return data[idx];
	}

	Vec2 operator+=(Vec2 v) {
		x += v.x;
		y += v.y;
		return *this;
	}
	Vec2 operator-=(Vec2 v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}
	Vec2 operator*=(Vec2 v) {
		x *= v.x;
		y *= v.y;
		return *this;
	}
	Vec2 operator/=(Vec2 v) {
		x /= v.x;
		y /= v.y;
		return *this;
	}

	Vec2 operator+=(float s) {
		x += s;
		y += s;
		return *this;
	}
	Vec2 operator-=(float s) {
		x -= s;
		y -= s;
		return *this;
	}
	Vec2 operator*=(float s) {
		x *= s;
		y *= s;
		return *this;
	}
	Vec2 operator/=(float s) {
		x /= s;
		y /= s;
		return *this;
	}

	Vec2 operator+(Vec2 v) const {
		return Vec2(x + v.x, y + v.y);
	}
	Vec2 operator-(Vec2 v) const {
		return Vec2(x - v.x, y - v.y);
	}
	Vec2 operator*(Vec2 v) const {
		return Vec2(x * v.x, y * v.y);
	}
	Vec2 operator/(Vec2 v) const {
		return Vec2(x / v.x, y / v.y);
	}

	Vec2 operator+(float s) const {
		return Vec2(x + s, y + s);
	}
	Vec2 operator-(float s) const {
		return Vec2(x - s, y - s);
	}
	Vec2 operator*(float s) const {
		return Vec2(x * s, y * s);
	}
	Vec2 operator/(float s) const {
		return Vec2(x / s, y / s);
	}

	bool operator==(Vec2 v) const {
		return x == v.x && y == v.y;
	}
	bool operator!=(Vec2 v) const {
		return x != v.x || y != v.y;
	}

	/// Absolute value
	Vec2 abs() const {
		return Vec2(std::abs(x), std::abs(y));
	}
	/// Negation
	Vec2 operator-() const {
		return Vec2(-x, -y);
	}
	/// Are all members real numbers?
	bool valid() const {
		return std::isfinite(x) && std::isfinite(y);
	}

	/// Modify vec to have unit length
	Vec2 normalize() {
		float n = norm();
		x /= n;
		y /= n;
		return *this;
	}
	/// Return unit length vec in the same direction
	Vec2 unit() const {
		float n = norm();
		return Vec2(x / n, y / n);
	}

	float norm_squared() const {
		return x * x + y * y;
	}
	float norm() const {
		return std::sqrt(norm_squared());
	}

	Vec2 range(float min, float max) const {
		if (!valid()) return Vec2();
		Vec2 r = *this;
		float range = max - min;
		while (r.x < min) r.x += range;
		while (r.x >= max) r.x -= range;
		while (r.y < min) r.y += range;
		while (r.y >= max) r.y -= range;
		return r;
	}

	union {
		struct {
			float x;
			float y;
		};
		float data[2] = {};
	};
};

inline Vec2 operator+(float s, Vec2 v) {
	return Vec2(v.x + s, v.y + s);
}
inline Vec2 operator-(float s, Vec2 v) {
	return Vec2(v.x - s, v.y - s);
}
inline Vec2 operator*(float s, Vec2 v) {
	return Vec2(v.x * s, v.y * s);
}
inline Vec2 operator/(float s, Vec2 v) {
	return Vec2(s / v.x, s / v.y);
}

/// Take minimum of each component
inline Vec2 hmin(Vec2 l, Vec2 r) {
	return Vec2(std::min(l.x, r.x), std::min(l.y, r.y));
}
/// Take maximum of each component
inline Vec2 hmax(Vec2 l, Vec2 r) {
	return Vec2(std::max(l.x, r.x), std::max(l.y, r.y));
}

/// 2D dot product
inline float dot(Vec2 l, Vec2 r) {
	return l.x * r.x + l.y * r.y;
}

inline std::string to_string(Vec2 const& v) {
	return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ")";
}

inline std::ostream& operator<<(std::ostream& out, Vec2 v) {
	out << "{" << v.x << "," << v.y << "}";
	return out;
}

// given line AB, CD,
// let result_bc = AB cross AC 
// let result_bd = AB cross AD 
// if result_bc * result_bd > 0 || (result_bc == 0 && result_bd == 0) 
//   AB not cross with CD
// else 
//   AB cross with CD 
inline bool is_line_cross(Vec2 pt_a, Vec2 pt_b, Vec2 pt_c, Vec2 pt_d) {
	auto cross = [](Vec2 u, Vec2 v) {
		return u.x * v.y - u.y * v.x;
	};
	auto ab = pt_b - pt_a, ac = pt_c - pt_a, ad = pt_d - pt_a,
		cd = pt_d - pt_c, ca = pt_a - pt_c, cb = pt_b - pt_c;
	auto result_abc = cross(ab, ac), result_abd = cross(ab, ad),
		result_cda = cross(cd, ca), result_cdb = cross(cd, cb);

	if (std::max(pt_a.x, pt_b.x) < std::min(pt_c.x, pt_d.x)
		|| std::max(pt_c.x, pt_d.x) < std::min(pt_a.x, pt_b.x)
		|| std::max(pt_a.y, pt_b.y) < std::min(pt_c.y, pt_d.y)
		|| std::max(pt_c.y, pt_d.y) < std::min(pt_a.y, pt_b.y)) {
		return false;
	}

	if (result_abc * result_abd <= 0 && result_cda * result_cdb <= 0) {
		return true;
	} else {
		return false;
	}
};