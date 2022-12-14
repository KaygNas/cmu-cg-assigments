
#pragma once

#include <memory>

#include "../lib/mathlib.h"

class Transform {
public:
	Transform() = default;
	Transform(Vec3 t, Vec3 euler, Vec3 s) : translation(t), rotation(Quat::euler(euler)), scale(s) {
	}
	Transform(Vec3 t, Quat q, Vec3 s) : translation(t), rotation(q), scale(s) {
	}

	std::weak_ptr<Transform> parent;
	Vec3 translation = Vec3{0.0f, 0.0f, 0.0f};
	Quat rotation = Quat::euler(Vec3{0.0f, 0.0f, 0.0f});
	Vec3 scale = Vec3{1.0f, 1.0f, 1.0f};

	Mat4 local_to_parent() const;
	Mat4 parent_to_local() const;
	Mat4 local_to_world() const;
	Mat4 world_to_local() const;
};

bool operator!=(const Transform& a, const Transform& b);
