#include "pipeline.h"

#include "framebuffer.h"
#include "sample_pattern.h"

#include "../lib/log.h"
#include "../lib/mathlib.h"
#include <list>
#include <iostream>
#include <set>
#include <algorithm>


template< PrimitiveType primitive_type, class Program, uint32_t flags >
void Pipeline< primitive_type, Program, flags >::run(
	std::vector< Vertex > const& vertices,
	typename Program::Parameters const& parameters,
	Framebuffer* framebuffer_) {
	//Framebuffer must be non-null:
	assert(framebuffer_);
	auto& framebuffer = *framebuffer_;

	//A1T7: sample loop
	//TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	// This will probably involve inserting a loop of the form:
	//     std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//     for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//  around some subset of the code.
	// You will also need to transform the input and output of the rasterize_* functions to
	//   account for the fact they deal with pixels centered at (0.5,0.5).

	std::vector< ShadedVertex > shaded_vertices;
	shaded_vertices.reserve(vertices.size());

	//--------------------------
	//shade vertices:
	for (auto const& v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex(parameters, v.attributes, &sv.clip_position, &sv.attributes);
		shaded_vertices.emplace_back(sv);
	}

	//--------------------------
	//assemble + clip + homogeneous divide vertices:
	std::vector< ClippedVertex > clipped_vertices;

	//reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		//clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		//clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	//helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const& sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	//actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line(shaded_vertices[i], shaded_vertices[i + 1], emit_vertex);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle(shaded_vertices[i], shaded_vertices[i + 1], shaded_vertices[i + 2], emit_vertex);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}


	//--------------------------
	//rasterize primitives:

	std::vector< Fragment > fragments;

	//helper used to put output of rasterization functions into fragments:
	auto emit_fragment = [&](Fragment const& f) {
		fragments.emplace_back(f);
	};
	//actually do rasterization:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
			rasterize_line(clipped_vertices[i], clipped_vertices[i + 1], emit_fragment);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
			rasterize_triangle(clipped_vertices[i], clipped_vertices[i + 1], clipped_vertices[i + 2], emit_fragment);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	//depth test + shade + blend fragments:
	uint32_t out_of_range = 0; //check if rasterization produced fragments outside framebuffer (indicates something is wrong with clipping)
	for (auto const& f : fragments) {

		//fragment location (in pixels):
		int32_t x = (int32_t)std::floor(f.fb_position.x);
		int32_t y = (int32_t)std::floor(f.fb_position.y);

		//if clipping is working properly, this condition shouldn't be needed;
		//however, it prevents crashes while you are working on your clipping functions,
		//so we suggest leaving it in place:
		if (x < 0 || (uint32_t)x >= framebuffer.width || y < 0 || (uint32_t)y >= framebuffer.height) {
			++out_of_range;
			continue;
		}

		//local names that refer to destination sample in framebuffer:
		float& fb_depth = framebuffer.depth_at(x, y, 0);
		Spectrum& fb_color = framebuffer.color_at(x, y, 0);

		//depth test:
		if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
			//"Always" means the depth test always passes.
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
			//"Never" means the depth test never passes.
			continue; //discard this fragment
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
			//"Less" means the depth test passes when the new fragment has depth less than the stored depth.
			//A1T4: Depth_Less
			//TODO: implement depth test!
		} else {
			static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
		}

		//if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
		if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
			fb_depth = f.fb_position.z;
		}

		//shade fragment:
		ShadedFragment sf;
		sf.fb_position = f.fb_position;
		Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

		//write color to framebuffer if color writes aren't disabled:
		if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {

			//blend fragment:
			if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
				fb_color = sf.color;
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
				//A1T4: Blend_Add
				//TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
				fb_color = sf.color; //<-- replace this line
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
				//A1T4: Blend_Over
				//TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
				fb_color = sf.color; //<-- replace this line
			} else {
				static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
			}
		}
	}

	if (out_of_range > 0) {
		if constexpr (primitive_type == PrimitiveType::Lines) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely wrong with the clip_line function.", out_of_range);
		} else if constexpr (primitive_type == PrimitiveType::Triangles) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely wrong with the clip_triangle function.", out_of_range);
		}
	}



}

//-------------------------------------------------------------------------
//clipping functions

//helper to interpolate between vertices:
template< PrimitiveType p, class P, uint32_t F >
auto Pipeline< p, P, F >::lerp(ShadedVertex const& a, ShadedVertex const& b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  va, vb: endpoints of line
 *  emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 *
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::clip_line(
	ShadedVertex const& va, ShadedVertex const& vb,
	std::function< void(ShadedVertex const&) > const& emit_vertex
) {
	//Determine portion of line over which:
	// pt = (b-a) * t + a
	// pt is a point on line
	// -pt.w <= pt.x <= pt.w
	// -pt.w <= pt.y <= pt.w
	// -pt.w <= pt.z <= pt.w

	//... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		//restrict range such that:
		//l + t * dl <= r + t * dr
		//re-arranging:
		// l - r <= t * (dr - dl)
		if (dr == dl) {
			//want: l - r <= 0
			if (l - r > 0.0f) {
				//works for none of range, so make range empty:
				min_t = 1.0f; max_t = 0.0f;
			}
		} else if (dr > dl) {
			//since dr - dl is positive:
			//want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { //dr < dl
			//since dr - dl is negative:
			//want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	//local names for clip positions and their difference:
	Vec4 const& a = va.clip_position;
	Vec4 const& b = vb.clip_position;
	Vec4 const ba = b - a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.x, ba.x);
	clip_range(a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.y, ba.y);
	clip_range(a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.z, ba.z);
	clip_range(a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va, vb, min_t);
			//don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) out.attributes = va.attributes;
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va, vb, max_t);
			//don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) out.attributes = va.attributes;
			emit_vertex(out);
		}
	}
}


/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  va, vb, vc: vertices of triangle
 *  emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 *
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::clip_triangle(
	ShadedVertex const& va, ShadedVertex const& vb, ShadedVertex const& vc,
	std::function< void(ShadedVertex const&) > const& emit_vertex
) {
	//A1T3: clip_triangle

	auto print_vec4 = [](auto prefix, Vec4 const& v) {
		std::cout << "\n" << prefix << ": " << v;
	};
	auto print_vertices = [&](std::vector<ShadedVertex>& vertices) {
		std::cout << "\nvertices:";
		int i = 0;
		for (ShadedVertex& vertex : vertices) {
			print_vec4(i, vertex.clip_position);
			++i;
		}
	};



	// ------------
	// clip triangle to produce more tirangle
	struct Plane {
		std::function<float(Vec4)> f;
		Vec4 n;
		Vec4 q;
	};
	float w = va.clip_position.w;
	float lw, rw, bw, tw, fw, nw;
	lw = bw = fw = 0 * w;
	rw = tw = nw = 1 * w;

	Vec4 lq(lw, 0.0f, 0.0f, w);
	Vec4 rq(rw, 0.0f, 0.0f, w);
	Vec4 bq(0.0f, bw, 0.0f, w);
	Vec4 tq(0.0f, tw, 0.0f, w);
	Vec4 fq(0.0f, 0.0f, fw, w);
	Vec4 nq(0.0f, 0.0f, nw, w);

	Vec4 ln(-lw, 0.0f, 0.0f, w);
	Vec4 rn(-rw, 0.0f, 0.0f, w);
	Vec4 bn(0.0f, -bw, 0.0f, w);
	Vec4 tn(0.0f, -tw, 0.0f, w);
	Vec4 fn(0.0f, 0.0f, -fw, w);
	Vec4 nn(0.0f, 0.0f, -nw, w);

	auto lf = [&](Vec4 pt) {
		return -pt.x + lw;
	};
	auto rf = [&](Vec4 pt) {
		return pt.x - rw;
	};
	auto bf = [&](Vec4 pt) {
		return -pt.y + bw;
	};
	auto tf = [&](Vec4 pt) {
		return pt.y - tw;
	};
	auto ff = [&](Vec4 pt) {
		return -pt.z + fw;
	};
	auto nf = [&](Vec4 pt) {
		return pt.z - nw;
	};

	std::vector <Plane> planes{
		 Plane{ lf, fn, fq },
		 Plane{ rf, rn, rq },
		 Plane{ bf, bn, bq },
		 Plane{ tf, tn, tq },
		 Plane{ ff, fn, fq },
		 Plane{ nf, nn, nq }
	};

	// for each of six planes do
	// if triangle entirely outside of plane then
	//   break
	// else if triangle spans plane then
	//   clip triangle
	//   if quadrilateral is left then
	//      break into two triangles
	std::vector<ShadedVertex> vertices{ va, vb, vc };

	for (Plane& plane : planes) {
		for (int i = 0; i + 2 < vertices.size();i += 3) {
			ShadedVertex va_ = vertices[i];
			ShadedVertex vb_ = vertices[i + 1];
			ShadedVertex vc_ = vertices[i + 2];
			auto fa = plane.f(va_.clip_position);
			auto fb = plane.f(vb_.clip_position);
			auto fc = plane.f(vc_.clip_position);
			// is entirely outside
			if (fa > 0 && fb > 0 && fc > 0) {
				std::cout << "\n------- vertices cleared ------\nvertices is entirely outside;";
				vertices.clear();
				std::cout << "\nvertices.size(): " << vertices.size() << "\n";
				break;
			}
			// is entirely inside
			else if (fa <= 0 && fb <= 0 && fc <= 0) {
			} else {

				if (fa * fc >= 0) {
					std::swap(fb, fc);
					std::swap(vb_, vc_);
					std::swap(fa, fb);
					std::swap(va_, vb_);
				} else if (fb * fc >= 0) {
					std::swap(fa, fc);
					std::swap(va_, vc_);
					std::swap(fa, fb);
					std::swap(va_, vb_);
				}

				auto intersection_vertex = [&](ShadedVertex va, ShadedVertex vb) {
					// t = n * a + D / (n * (a - b)); D = -n * q;
					auto D = dot(-plane.n, plane.q);
					auto t = (dot(plane.n, va.clip_position) + D) / dot(plane.n, (va.clip_position - vb.clip_position));
					return lerp(va, vb, t);
				};
				auto A = intersection_vertex(va_, vc_);
				auto B = intersection_vertex(vb_, vc_);

				if (fc <= 0) {
					vertices[i] = A;
					vertices[i + 1] = B;
					vertices[i + 2] = vc_;
				} else {
					vertices[i] = va_;
					vertices[i + 1] = vb_;
					vertices[i + 2] = A;
					vertices.insert(vertices.begin() + i + 3, { vb_,B,A });
				}

				std::cout << "\n cutting: " << i / 3 << " q:" << plane.q;
				std::cout << "\na: " << va_.clip_position << " b: " << vb_.clip_position << " c: " << vc_.clip_position;
				std::cout << "\nfa: " << fa << " fb: " << fb << " fc: " << fc;
				std::cout << "\nA: " << A.clip_position << " B: " << B.clip_position;

				print_vertices(vertices);
			}
		}
	}

	// -------------
	// compact triangles and minimize the number of triangles
	auto join_polygon = [&](std::vector<ShadedVertex>& a, std::vector<ShadedVertex>& b) {
		std::cout << "\n\na ----";
		print_vertices(a);
		std::cout << "\nb ----";
		print_vertices(b);

		std::vector<ShadedVertex> vertices;
		for (int i = 1; i <= a.size(); ++i) {
			if (vertices.size() > 0) break;

			ShadedVertex va_1 = a[i - 1];
			ShadedVertex va_2 = a[i % a.size()];
			for (int j = 1; j <= b.size(); ++j) {
				ShadedVertex vb_1 = b[j - 1];
				ShadedVertex vb_2 = b[j % b.size()];
				if (va_1.clip_position == vb_2.clip_position && va_2.clip_position == vb_1.clip_position) {
					print_vec4("va_1", va_1.clip_position);
					print_vec4("va_2", va_2.clip_position);
					print_vec4("vb_1", vb_1.clip_position);
					print_vec4("vb_2", vb_2.clip_position);
					auto emplace_vertex = [&](ShadedVertex& vertex) {
						if (vertices.size() < 2) {
							vertices.emplace_back(vertex);
							return;
						}
						// check is the vertex is on the same line
						// if is, replace the vertex
						// else emplace
						Vec4& a = vertices[vertices.size() - 2].clip_position;
						Vec4& b = vertices[vertices.size() - 1].clip_position;
						Vec4& c = vertex.clip_position;
						if ((b - a).unit() == (c - a).unit()) {
							print_vec4("replacing", b);
							print_vec4("to", c);
							vertices[vertices.size() - 1] = vertex;
						} else {
							vertices.emplace_back(vertex);
						}
					};
					for (int n = 0; n < i; ++n) {
						print_vec4("a[n]", a[n].clip_position);
						emplace_vertex(a[n]);
					}
					for (int k = (j + 1) % b.size(); k < (j + 1 >= b.size() ? j - 1 : b.size()); ++k) {
						print_vec4("b[k]", b[k].clip_position);
						emplace_vertex(b[k]);
					}
					for (int l = i; l < a.size(); ++l) {
						print_vec4("a[l]", a[l].clip_position);
						emplace_vertex(a[l]);
					}

					break;
				}
			}
		}
		return vertices;
	};

	auto slice_vector = [](std::vector<ShadedVertex>& vector, int start, int end) {
		if (vector.size() == 0) {
			std::vector<ShadedVertex> new_vector;
			return new_vector;
		}

		start = std::clamp(start, 0, (int)vector.size() - 1);
		end = std::clamp(end, start, (int)vector.size());
		std::cout << "\n start: " << start << " end: " << end;
		std::vector<ShadedVertex> new_vector(vector.begin() + start, vector.begin() + end);
		return new_vector;
	};
	std::vector<ShadedVertex> polygon = slice_vector(vertices, 0, 3);
	std::vector<ShadedVertex> vertices_candidates = slice_vector(vertices, 3, INFINITY);
	int last_size = INFINITY;
	while (vertices_candidates.size() != last_size) {
		last_size = vertices_candidates.size();
		for (int i = 2; i < vertices_candidates.size();i += 3) {
			std::vector<ShadedVertex> next_polygon(vertices_candidates.begin() + i - 2, vertices_candidates.begin() + i + 1);
			std::vector<ShadedVertex> new_polygon = join_polygon(polygon, next_polygon);
			if (new_polygon.size() > 0) {
				polygon = new_polygon;
				vertices_candidates.erase(vertices_candidates.begin() + i - 2, vertices_candidates.begin() + i + 1);
				std::cout << "\nvertices_candidates.size() " << vertices_candidates.size();
				std::cout << "\nlast_size " << last_size;
				print_vertices(vertices_candidates);
				break;
			}
		}
	}

	// split polygon into mutilple triangles
	std::vector<ShadedVertex> triangles;
	if (polygon.size() > 0) {
		ShadedVertex v_origin = polygon.front();
		for (int i = 2; i < polygon.size(); ++i) {
			triangles.emplace_back(v_origin);
			triangles.emplace_back(polygon[i - 1]);
			triangles.emplace_back(polygon[i]);
		}
	}

	for (ShadedVertex& vertex : triangles) {
		emit_vertex(vertex);
	}
}


//-------------------------------------------------------------------------
//rasterization functions


/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *   a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 *       (x+0.5,y+1)
 *       /        \
 *   (x,y+0.5)  (x+1,y+0.5)
 *       \        /
 *        (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points.
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 *
 */

template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::rasterize_line(
	ClippedVertex const& va, ClippedVertex const& vb,
	std::function< void(Fragment const&) > const& emit_fragment
) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	//A1T2: rasterize_line
	// pt = (vb - va) * t + va
	// pt.x = (vb.x - va.x) * t + va.x
	// pt.y = (vb.y - va.y) * t + va.y
	// slop_x = vb.x - va.x, slop_y = vb.y - va.y
	// if |slop_x| >= |slop_y| then:
	//   increment at x
	// else:
	//   increment at y

	auto a = va.fb_position.xy();
	auto b = vb.fb_position.xy();
	float slop_x = b.x - a.x;
	float slop_y = b.y - a.y;

	auto fy = [&](float x) {
		float t = (x - a.x) / slop_x;
		return slop_y * t + a.y;
	};
	auto fx = [&](float y) {
		float t = (y - a.y) / slop_y;
		return slop_x * t + a.x;
	};
	auto round_xy = [](float v) {
		float fv = std::floor(v);
		float cv = std::ceil(v);
		float offset = v > 0 ? 0.5 : -0.5;
		if (fv <= v && v < cv) {
			return fv + offset;
		} else {
			return cv + offset;
		}
	};
	auto iterateRange = [](
		float const& start, float const& end, float const& step_,
		std::function< void(float const&) > const& emit_value
		) {
			assert(step_ > 0 && "step should greater than 0");

			bool const incremental = end - start > 0;
			float const step = incremental ? step_ : -step_;
			float value = start + step;
			while (incremental ? value < end : value > end) {
				emit_value(value);
				value = value + step;
			}
	};
	auto make_fragment = [&](float const x, float const y, float const z) {
		Fragment mid;
		mid.fb_position = Vec3(x, y, z);
		mid.attributes = va.attributes;
		mid.derivatives.fill(Vec2(0.0f, 0.0f));
		return mid;
	};
	auto interpolate_z = [&va, &vb](float const x_i) {
		return (x_i - va.fb_position.x) / (vb.fb_position.x - va.fb_position.x) * (vb.fb_position.z - va.fb_position.z) + va.fb_position.z;
	};
	float r_ax = round_xy(a.x), r_bx = round_xy(b.x);
	float r_ay = round_xy(a.y), r_by = round_xy(b.y);
	auto is_a_enterexit = is_enterexit_diamond(a, b, Vec2(std::floor(a.x), std::floor(a.y)));
	auto is_b_enterexit = is_enterexit_diamond(a, b, Vec2(std::floor(b.x), std::floor(b.y)));
	if (is_a_enterexit) {
		emit_fragment(make_fragment(r_ax, r_ay, va.fb_position.z));
	}
	if (is_b_enterexit) {
		emit_fragment(make_fragment(r_bx, r_by, vb.fb_position.z));
	}
	if (abs(slop_x) >= abs(slop_y)) {
		assert(slop_x != 0.0f && "slop_x should not be 0");
		iterateRange(
			r_ax, r_bx, 1,
			[&](float x_i) {
				auto y_i = fy(x_i);
				emit_fragment(make_fragment(x_i, round_xy(y_i), interpolate_z(x_i)));
			});
	} else {
		assert(slop_y != 0.0f && "slop_y should not be 0");
		iterateRange(
			r_ay, r_by, 1,
			[&](float y_i) {
				auto x_i = fx(y_i);
				emit_fragment(make_fragment(round_xy(x_i), y_i, interpolate_z(x_i)));
			});
	}

}


/*
 *
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  (x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Screen: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 *  The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/dfb_position.x attributes[i].x
 *    derivatives[i].y = d/dfb_position.y attributes[i].y
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::rasterize_triangle(
	ClippedVertex const& va, ClippedVertex const& vb, ClippedVertex const& vc,
	std::function< void(Fragment const&) > const& emit_fragment
) {
	//NOTE: it is okay to restructure this function to allow these tasks to use the
	// same code paths. Be aware, however, that all of them need to remain working!
	// (e.g., if you break Flat while implementing Correct, you won't get points
	//  for Flat.)
	if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
		//A1T3: flat triangles
		//TODO: rasterize triangle (see block comment above this function).

		//As a placeholder, here's code that draws some lines:
		//(remove this and replace it with a real solution)
		Pipeline< PrimitiveType::Lines, P, flags >::rasterize_line(va, vb, emit_fragment);
		Pipeline< PrimitiveType::Lines, P, flags >::rasterize_line(vb, vc, emit_fragment);
		Pipeline< PrimitiveType::Lines, P, flags >::rasterize_line(vc, va, emit_fragment);
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Screen) {
		//A1T5: screen-space smooth triangles
		//TODO: rasterize triangle (see block comment above this function).

		//As a placeholder, here's code that calls the Flat interpolation version of the function:
		//(remove this and replace it with a real solution)
		Pipeline< PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Flat >::rasterize_triangle(va, vb, vc, emit_fragment);
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		//A1T5: perspective correct triangles
		//TODO: rasterize triangle (block comment above this function).

		//As a placeholder, here's code that calls the Screen-space interpolation function:
		//(remove this and replace it with a real solution)
		Pipeline< PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Screen >::rasterize_triangle(va, vb, vc, emit_fragment);
	}
}


//-------------------------------------------------------------------------
//compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Screen >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct >;
