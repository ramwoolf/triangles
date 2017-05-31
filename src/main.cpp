#include <iostream>
#include <cmath>
#include <array>
#include <limits>
#include <algorithm>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>

#include "/usr/include/tinyxml2.h"

namespace VolkovTest 
{
	constexpr int input_array_dimension = 9;
	constexpr int input_chunk_dimension = 3;

	enum class POINT2SEGMENT_RELATION {LEFT, RIGHT};

	struct Point3D;
	using Vector3D = Point3D;
	using ScalarTriplet = Point3D;
	using LineParameterSegment = std::pair<double, double>;
	using Segment = std::pair<Point3D, Point3D>;

	double dot_product(const Vector3D&, const Vector3D&);
	double dot_product(const Point3D&, const Point3D&, const Point3D&, const Point3D&);
	Vector3D cross_product(const Vector3D&, const Vector3D&);
	Vector3D cross_product(const Point3D&, const Point3D&, const Point3D&, const Point3D&);

	using lim_double = std::numeric_limits<double>;

    bool is_two_doubles_equal(double x, double y, int ulp = 6)
	{
    	return std::abs(x - y) < lim_double::epsilon() * std::abs(x + y) * ulp
           || std::abs(x - y) < lim_double::min();
	}

	struct Point3D
	{
		double x;
		double y;
		double z;

		Point3D(double _x = 0.0, double _y = 0.0, double _z = 0.0)
		: x(_x), y(_y), z(_z)
		{}

		Point3D(Point3D&&) = default;
		Point3D& operator=(Point3D&&) = default;

		Point3D(const Point3D&) = default;
		Point3D& operator=(const Point3D&) = default;

		bool is_null() const
		{
			return is_two_doubles_equal(x, 0.0) 
				&& is_two_doubles_equal(y, 0.0)
				&& is_two_doubles_equal(z, 0.0);
		}

		bool all_have_same_sign() const
		{
			if (is_null()) return false;
			if ((x < 0.0 && y < 0.0 && z < 0.0) || (x > 0.0 && y > 0.0 && z > 0.0)) return true;
			return false;
		}

		Point3D& operator+=(const Point3D& rhs)
		{
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
			return *this;
		}

		Point3D& operator-=(const Point3D& rhs)
		{
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
			return *this;
		}

		double& operator[](int index)
		{
			switch (index)
			{
				case 0: return x;
				case 1: return y;
				case 2:	return z;
				default: throw std::range_error("Bad index for 3D element");
			}
		}

		double operator[](int index) const
		{
			switch (index)
			{
				case 0: return x;
				case 1: return y;
				case 2:	return z;
				default: throw std::range_error("Bad index for 3D element");
			}
		}

		friend Point3D operator+(Point3D&& lhs, const Point3D& rhs)
		{
			lhs += rhs;
			return std::move(lhs);
		}

		friend Point3D operator-(Point3D&& lhs, const Point3D& rhs)
		{
			lhs -= rhs;
			return std::move(lhs);
		}

		friend Point3D operator-(const Point3D& lhs, const Point3D& rhs)
		{
			Point3D result;
			for (int i = 0; i < input_chunk_dimension; ++i)
			{
				result[i] = lhs[i] - rhs[i];
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& output, const Point3D& arg)
		{
			std::cout << "Point3D" << std::endl;
			std::cout << arg.x << ", " << arg.y << ", " << arg.z << std::endl;
			return output;
		}
	};

	struct Triangle3D
	{
		Triangle3D(Point3D _a, Point3D _b, Point3D _c)
		: a(_a), b(_b), c(_c)
		{}

		Point3D a;
		Point3D b;
		Point3D c;

		Point3D& operator[](int index)
		{
			switch (index)
			{
				case 0: return a;
				case 1:	return b;
				case 2:	return c;
				default: throw std::range_error("Bad index for triangle vertex");
			}
		}

		Point3D operator[](int index) const
		{
			switch (index)
			{
				case 0:	return a;
				case 1:	return b;
				case 2:	return c;
				default: throw std::range_error("Bad index for triangle vertex");
			}
		}
	};


	// Plane represented by equation: normal_vector * point + d = 0
	class Plane3D
	{
		public:
			Plane3D(const Triangle3D& triangle)
			{
				normal_vector = cross_product(triangle.a, triangle.b, triangle.a, triangle.c);
				d = -normal_vector.x * triangle.b.x - normal_vector.y * triangle.b.y - normal_vector.z * triangle.b.z;
			}

			const Vector3D& get_normal_vector() const
			{
				return normal_vector;
			}

			const double get_plane_shift() const
			{
				return d;
			}

		private:
			Vector3D normal_vector;
			double d;
	};

	struct TwoElemsAndOpposite
	{
		std::pair<double, double> two_similar_params;
		double opposite_param;
	};

	double dot_product(const Vector3D& lhs, const Vector3D& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
	}

	double dot_product(const Point3D& lhs_start, const Point3D& lhs_end, const Point3D& rhs_start, const Point3D& rhs_end)
	{
		return (lhs_end.x - lhs_start.x) * (rhs_end.x - rhs_start.x) 
			 + (lhs_end.y - lhs_start.y) * (rhs_end.y - rhs_start.y)  
			 + (lhs_end.z - lhs_start.z) * (rhs_end.z - rhs_start.z);
	}

	Vector3D cross_product(const Vector3D& lhs, const Vector3D& rhs)
	{
		return Vector3D(lhs.y * rhs.z - lhs.z * rhs.y, lhs.x * rhs.z - lhs.z * rhs.x, lhs.y * rhs.x - lhs.x * rhs.y);
	}

	Vector3D cross_product(const Point3D& lhs_start, const Point3D& lhs_end, const Point3D& rhs_start, const Point3D& rhs_end)
	{
		return Vector3D( (lhs_end.y - lhs_start.y) * (rhs_end.z - rhs_start.z) - (lhs_end.z - lhs_start.z) * (rhs_end.y - rhs_start.y),
						 (lhs_end.x - lhs_start.x) * (rhs_end.z - rhs_start.z) - (lhs_end.z - lhs_start.z) * (rhs_end.x - rhs_start.x), 
						 (lhs_end.y - lhs_start.y) * (rhs_end.x - rhs_start.x) - (lhs_end.x - lhs_start.x) * (rhs_end.y - rhs_start.y));
	}

	double triple_product(const Vector3D& a, const Vector3D& b, const Vector3D& c)
	{
		return dot_product(a, cross_product(b, c));
	}

	auto convert_plain_array_to_triplets(const double arg[])
	{
		return std::array<Point3D, input_chunk_dimension> {Point3D(arg[0], arg[1], arg[2]), Point3D(arg[3], arg[4], arg[5]), Point3D(arg[6], arg[7], arg[8])};
	}

	bool can_three_points_be_triangle(const auto& arg)
	{
		Vector3D first_vector{ arg[1].x - arg[0].x, arg[1].y - arg[0].y, arg[1].z - arg[0].z};
		Vector3D second_vector{arg[2].x - arg[0].x, arg[2].y - arg[0].y, arg[2].z - arg[0].z};
		if (first_vector.is_null() || second_vector.is_null()) return false;

		return !cross_product(arg[0], arg[1], arg[0], arg[2]).is_null();
	}

	ScalarTriplet get_signed_distances(const Triangle3D& triangle, const Plane3D& plane)
	{
		auto normal_vector = plane.get_normal_vector();
		auto d = plane.get_plane_shift();

		ScalarTriplet result;
		for (int i = 0; i < input_chunk_dimension; ++i)
		{
			result[i] = dot_product(normal_vector, triangle[i]) + d;
		}
		return result;
	}

	int get_point_coord_index(Vector3D normal_vector)
	{
		for (int i = 0; i < input_chunk_dimension; ++i)
		{
			normal_vector[i] = std::abs(normal_vector[i]);
		}
		auto max_normal_component = std::max({normal_vector.x, normal_vector.y, normal_vector.z});

		if 		(max_normal_component == normal_vector.x) return 0;
		else if (max_normal_component == normal_vector.y) return 1;
		else if (max_normal_component == normal_vector.z) return 2;
		else throw std::range_error("Wrong index of max element"); // impossible
	}

	double calc_line_param(	const double p_from_twice,
							const double p_alone,
							const double d_from_twice,
							const double d_alone)
	{
		return p_from_twice + (d_from_twice * (p_alone - p_from_twice)) / (d_from_twice - d_alone);
	}

	TwoElemsAndOpposite sort_by_sign(const ScalarTriplet& criteria, const ScalarTriplet& values)
	{
		TwoElemsAndOpposite result;

		auto criteria_x_is_null = is_two_doubles_equal(criteria.x, 0.0);
		auto criteria_y_is_null = is_two_doubles_equal(criteria.y, 0.0);
		auto criteria_z_is_null = is_two_doubles_equal(criteria.z, 0.0);

		if (   (criteria.x > 0 && criteria.y > 0 && criteria.z < 0) 
			|| (criteria.x < 0 && criteria.y < 0 && criteria.z > 0)
			|| (criteria_x_is_null && criteria_y_is_null && !criteria_z_is_null)
			|| (!criteria_x_is_null && !criteria_y_is_null && criteria_z_is_null))
		{
			result.two_similar_params.first = values.x;
			result.two_similar_params.second = values.y;
			result.opposite_param = values.z;
		}
		else if (  (criteria.x > 0 && criteria.y < 0 && criteria.z > 0) 
				|| (criteria.x < 0 && criteria.y > 0 && criteria.z < 0)
				|| (criteria_x_is_null && !criteria_y_is_null && criteria_z_is_null)
				|| (!criteria_x_is_null && criteria_y_is_null && !criteria_z_is_null))
		{
			result.two_similar_params.first = values.x;
			result.two_similar_params.second = values.z;
			result.opposite_param = values.y;
		}
		else if (  (criteria.x > 0 && criteria.y < 0 && criteria.z < 0) 
				|| (criteria.x < 0 && criteria.y > 0 && criteria.z > 0)
				|| (criteria_x_is_null && !criteria_y_is_null && !criteria_z_is_null)
				|| (!criteria_x_is_null && criteria_y_is_null && criteria_z_is_null))
		{
			result.two_similar_params.first = values.y;
			result.two_similar_params.second = values.z;
			result.opposite_param = values.x;
		}
		else
		{
			throw std::logic_error("Wrong signs set");
		}
		return result;
	}

	LineParameterSegment get_line_param_segment(const ScalarTriplet& signed_distances, const ScalarTriplet& intersect_repr)
	{
		LineParameterSegment result;

		auto repr = sort_by_sign(signed_distances, intersect_repr);
		auto dist = sort_by_sign(signed_distances, signed_distances);

		result.first  = calc_line_param(repr.two_similar_params.first, repr.opposite_param, 
										dist.two_similar_params.first, dist.opposite_param);
		result.second = calc_line_param(repr.two_similar_params.second, repr.opposite_param,
										dist.two_similar_params.second, dist.opposite_param);
		if (result.first > result.second) std::swap(result.first, result.second);
		return result;
	}

	bool are_segments_intersected(const LineParameterSegment& first_segment, const LineParameterSegment& second_segment)
	{
		return (first_segment.first >= second_segment.first && first_segment.first <= second_segment.second)
				|| (first_segment.second >= second_segment.first && first_segment.second <= second_segment.second);
	}

	bool check_isect_in_segment_bounds(const Point3D& point, const Segment& segment)
	{
		auto main_vector = segment.first - segment.second;
		auto cross_test_vector = segment.first - point;
		auto first_test_vector = point - segment.first;
		auto second_test_vector = point - segment.second;

		return cross_product(main_vector, cross_test_vector).is_null() && dot_product(first_test_vector, second_test_vector) <= 0;
	}

	bool check_segments_intersection(const Segment& a, const Segment& b)
	{
		auto result = false;
		Vector3D fictive_vector(1, 1, 1);
		auto main_vector = a.second - a.first;
		auto first_test_vector = b.first - a.first;
		auto second_test_vector = b.second - a.first;

		if (  triple_product(fictive_vector, main_vector, first_test_vector) 
			* triple_product(fictive_vector, main_vector, second_test_vector) < 0)
		{
			result |= true;
		}

		if (result) return result;

		result = result || check_isect_in_segment_bounds(a.first, b);
		result = result || check_isect_in_segment_bounds(a.second, b);
		result = result || check_isect_in_segment_bounds(b.first, a);
		result = result || check_isect_in_segment_bounds(b.second, a);

		return result;
	}

	POINT2SEGMENT_RELATION detect_relation(const Point3D& point, const Point3D& start, const Point3D& end)
	{
		Vector3D fictive_vector(1, 1, 1);
		return triple_product(fictive_vector, std::move(end) - start, std::move(point) - start) > 0 ?
		 	POINT2SEGMENT_RELATION::LEFT :
			POINT2SEGMENT_RELATION::RIGHT;
	}

	bool point_in_triangle_test(const Point3D& point, const Triangle3D& triangle)
	{
		auto relation_ab = detect_relation(point, triangle.a, triangle.b);
		auto relation_bc = detect_relation(point, triangle.b, triangle.c);
		auto relation_ca = detect_relation(point, triangle.c, triangle.a);
		return relation_ab == relation_bc && relation_bc == relation_ca;
	}

	bool co_planar_triangles_test_intersection(const Triangle3D& first_triangle, const Triangle3D& second_triangle, const Plane3D& plane)
	{
		auto result = false;
		for (int i = 0; i < input_chunk_dimension; ++i)
		{
			auto next_i = (i + 1) % input_chunk_dimension;
			for (int j = 0; j < input_chunk_dimension; ++j)
			{
				auto next_j = (j + 1) % input_chunk_dimension;
				Segment a = std::make_pair(first_triangle[i], first_triangle[next_i]);
				Segment b = std::make_pair(second_triangle[j], second_triangle[next_j]);
				result = result || check_segments_intersection(a, b);
			}
		}

		if (result) return result;

		for (int i = 0; i < input_chunk_dimension; i++)
		{
			result = result || point_in_triangle_test(first_triangle[i], second_triangle);
			result = result || point_in_triangle_test(second_triangle[i], first_triangle);
		}

		return result;
	}

	bool test_line_triangles_intersection(const Triangle3D& first_triangle, const Triangle3D& second_triangle,
											const Plane3D& first_plane, const Plane3D& second_plane,
											const ScalarTriplet& signed_distances_t1_p2, const ScalarTriplet& signed_distances_t2_p1)
	{
		auto intersect_line_dir_vector = cross_product(first_plane.get_normal_vector(), second_plane.get_normal_vector());
		auto point_coord_index = get_point_coord_index(intersect_line_dir_vector);	

		ScalarTriplet scalar_intersect_repr_t1, scalar_intersect_repr_t2;
		for (int i = 0; i < input_chunk_dimension; ++i)
		{
			scalar_intersect_repr_t1[i] = first_triangle[i][point_coord_index];
			scalar_intersect_repr_t2[i] = second_triangle[i][point_coord_index];
		}

		try
		{
			LineParameterSegment first_segment = get_line_param_segment(signed_distances_t1_p2, scalar_intersect_repr_t1);
			LineParameterSegment second_segment = get_line_param_segment(signed_distances_t2_p1, scalar_intersect_repr_t2);

			return are_segments_intersected(first_segment, second_segment);
		}
		catch(...)
		{
			return false;
		}
	}

	bool check_triangles_intersection(const double a[], const double b[])
	{
		auto first_coord_triplet = convert_plain_array_to_triplets(a);
		auto second_coord_triplet = convert_plain_array_to_triplets(b);

		if (!can_three_points_be_triangle(first_coord_triplet))
		{
			return false; // first three points cannot be triangle
		}

		if (!can_three_points_be_triangle(second_coord_triplet))
		{
			return false; // second three points cannot be triangle
		}

		Triangle3D first_triangle(first_coord_triplet[0], first_coord_triplet[1], first_coord_triplet[2]);
		Triangle3D second_triangle(second_coord_triplet[0], second_coord_triplet[1], second_coord_triplet[2]);

		Plane3D first_plane(first_triangle);
		Plane3D second_plane(second_triangle);

		auto signed_distances_t1_p2 = get_signed_distances(first_triangle, second_plane);
		auto signed_distances_t2_p1 = get_signed_distances(second_triangle, first_plane);

		if (signed_distances_t1_p2.is_null())
		{
			// triangles are co-planar
			return co_planar_triangles_test_intersection(first_triangle, second_triangle, first_plane);
		}
		else
		{
			// triangles aren't co-planar
			if (signed_distances_t1_p2.all_have_same_sign() || signed_distances_t2_p1.all_have_same_sign())
			{
				return false;
			}
			else
			{
				return test_line_triangles_intersection(first_triangle, second_triangle,
														 first_plane, second_plane,
														 signed_distances_t1_p2, signed_distances_t2_p1);
			}
		}
	}
}

class TestCase
{

}

int main(int argc, char const *argv[]) 
{
	XMLDocument xmlDoc;
	xmlDoc.LoadFile("TriTestBase.xml");

	XMLNode *suiteNode = xmlDoc.FirstChild();




	return 0;
}