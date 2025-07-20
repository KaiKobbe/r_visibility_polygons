
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <fmt/core.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>

using Kernel = CGAL::Epeck;
using Point = CGAL::Point_2<Kernel>;
using Polygon2WithHoles = CGAL::Polygon_with_holes_2<Kernel>;
using Polygon2 = CGAL::Polygon_2<Kernel>;

using Segment2 = Kernel::Segment_2;
using Traits_2 = CGAL::Arr_segment_traits_2<Kernel>;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2>;
using Halfedge_const_handle = Arrangement_2::Halfedge_const_handle;
using Face_handle = Arrangement_2::Face_handle;
using PointLocation = CGAL::Arr_naive_point_location<Arrangement_2>;

using Arr_face_extended_dcel = CGAL::Arr_face_extended_dcel<Traits_2, std::set<int>>;
using Ex_arrangement = CGAL::Arrangement_2<Traits_2, Arr_face_extended_dcel>;
using Arr_point_location = CGAL::Arr_trapezoid_ric_point_location<Arrangement_2>;

using Point_2 = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;

class Rectilinear_Visibility {
public:

    Rectilinear_Visibility(const Arrangement_2& arrangement)
    // Computes the r-visibility-polygon (without "antennas") for vertex guards in orthogonal polygons.
            : env(arrangement) {}

    Polygon2 compute_visibility(const Point_2& query_point, Halfedge_const_handle he) {

        // Decide which quadrants are seen by the query point
        auto pred = he->source();
        auto succ = he->next()->target();

        bool is_big_cone = CGAL::right_turn(pred->point(), query_point, succ->point());
        bool is_top_left = false;
        bool is_top_right = false;
        bool is_bottom_left = false;
        bool is_bottom_right = false;

        // Find affected quadrants
        if (pred->point().x() > query_point.x()) {
            if(is_big_cone) {is_bottom_left = true;} else {is_bottom_right = true;}
        }
        else if (pred->point().x() < query_point.x()) {
            if(is_big_cone) {is_top_right = true;} else {is_top_left = true;}
        }
        else if (pred->point().y() > query_point.y()) {
            if(is_big_cone) {is_bottom_right = true;} else {is_top_right = true;}
        }
        else if (pred->point().y() < query_point.y()) {
            if(is_big_cone) {is_top_left = true;} else {is_bottom_left = true;}
        }

        // Find canonical point set for each quadrant and compute chain out of them
        std::queue<Point_2> q1_points, q2_points, q3_points, q4_points;
        std::vector<Point_2> q1_chain, q2_chain, q3_chain, q4_chain;
        std::set<Halfedge_const_handle> processed_halfedges;

        for (auto he = env.halfedges_begin(); he != env.halfedges_end(); ++he) {

            if (processed_halfedges.find(he) != processed_halfedges.end()) {continue;}
            processed_halfedges.insert(he->twin());

            Segment_2 segment = he->curve();
            Point_2 source = segment.source();
            Point_2 target = segment.target();

            bool source_in_q1 = (source.x() > query_point.x() && source.y() > query_point.y());
            bool target_in_q1 = (target.x() > query_point.x() && target.y() > query_point.y());
            bool source_in_q2 = (source.x() < query_point.x() && source.y() > query_point.y());
            bool target_in_q2 = (target.x() < query_point.x() && target.y() > query_point.y());
            bool source_in_q3 = (source.x() < query_point.x() && source.y() < query_point.y());
            bool target_in_q3 = (target.x() < query_point.x() && target.y() < query_point.y());
            bool source_in_q4 = (source.x() > query_point.x() && source.y() < query_point.y());
            bool target_in_q4 = (target.x() > query_point.x() && target.y() < query_point.y());

            if ((is_big_cone && !is_bottom_left) || is_top_right) {
                process_segment(q1_points, source, target, query_point, source_in_q1, target_in_q1, is_vertical(segment), true, true);
            }
            if ((is_big_cone && !is_bottom_right) || is_top_left) {
                process_segment(q2_points, source, target, query_point, source_in_q2, target_in_q2, is_vertical(segment), false, true);
            }
            if ((is_big_cone && !is_top_right) || is_bottom_left) {
                process_segment(q3_points, source, target, query_point, source_in_q3, target_in_q3, is_vertical(segment), false, false);
            }
            if ((is_big_cone && !is_top_left) || is_bottom_right) {
                process_segment(q4_points, source, target, query_point, source_in_q4, target_in_q4, is_vertical(segment), true, false);
            }
        }

        sort_and_filter_points(q1_points, true, false);
        quadrant_chain(query_point, q1_chain, q1_points);
        std::reverse(q1_chain.begin(), q1_chain.end());

        sort_and_filter_points(q2_points, false, false);
        quadrant_chain(query_point, q2_chain, q2_points);

        sort_and_filter_points(q3_points, false, true);
        quadrant_chain(query_point, q3_chain, q3_points);
        std::reverse(q3_chain.begin(), q3_chain.end());

        sort_and_filter_points(q4_points, true, true);
        quadrant_chain(query_point, q4_chain, q4_points);

        // concatenate respective chains
        std::vector<Point_2> polygon = build_polygon_from_chains(query_point, q1_chain, q2_chain, q3_chain, q4_chain, is_big_cone, is_top_left, is_top_right, is_bottom_left, is_bottom_right);

        remove_duplicates(polygon);

        return make_polygon(polygon);
    }

private:

    const Arrangement_2& env;

    Polygon2 make_polygon(const std::vector<Point_2>& polygon) {
        Polygon2 final_polygon;
        for (const auto& point : polygon) {
            final_polygon.push_back(point);
        }
        return final_polygon;
    }

    void remove_duplicates(std::vector<Point_2>& points) {
        if (points.empty()) return;

        auto it = points.begin();
        while (it != points.end() - 1) {
            if (*it == *(it + 1)) {
                it = points.erase(it);
            } else {
                ++it;
            }
        }
    }

    std::vector<Point_2> build_polygon_from_chains(Point_2 query_point, std::vector<Point_2>& q1_chain, std::vector<Point_2>& q2_chain, std::vector<Point_2>& q3_chain, std::vector<Point_2>& q4_chain, bool is_big_cone, bool is_top_left, bool is_top_right, bool is_bottom_left, bool is_bottom_right) {
        if (is_big_cone && is_top_left) {
            // top_right -> top_left -> bottom_left -> qp
            q3_chain.push_back(query_point);
            return concatenate_vectors(q1_chain, q2_chain, q3_chain);

        } else if (is_big_cone && is_top_right) {
            // bottom_right -> top_right -> top_left -> qp
            q2_chain.push_back(query_point);
            return concatenate_vectors(q4_chain, q1_chain, q2_chain);

        } else if (is_big_cone && is_bottom_left) {
            // top_left -> bottom_left -> bottom_right -> qp
            q4_chain.push_back(query_point);
            return concatenate_vectors(q2_chain, q3_chain, q4_chain);

        } else if (is_big_cone && is_bottom_right) {
            // bottom_left -> bottom_right -> top_right -> qp
            q1_chain.push_back(query_point);
            return concatenate_vectors(q3_chain, q4_chain, q1_chain);
        }

        if (!is_big_cone && is_top_left) {
            q2_chain.push_back(query_point);
            return q2_chain;
        } else if (!is_big_cone && is_top_right) {
            q1_chain.push_back(query_point);
            return q1_chain;
        } else if (!is_big_cone && is_bottom_left) {
            q3_chain.push_back(query_point);
            return q3_chain;
        } else if (!is_big_cone && is_bottom_right) {
            q4_chain.push_back(query_point);
            return q4_chain;
        }
    }

    std::vector<Point_2> concatenate_vectors(const std::vector<Point_2>& vec1, const std::vector<Point_2>& vec2, const std::vector<Point_2>& vec3) {
        std::vector<Point_2> result;
        result.reserve(vec1.size() + vec2.size() + vec3.size());
        result.insert(result.end(), vec1.begin(), vec1.end());
        result.insert(result.end(), vec2.begin(), vec2.end());
        result.insert(result.end(), vec3.begin(), vec3.end());

        return result;
    }

    void quadrant_chain(Point_2 query_point, std::vector<Point_2>& q_chain, std::queue<Point_2>& q_points) {
        if (q_points.empty()) {return;}
        q_chain.push_back(q_points.front());
        q_points.pop();
        while (!q_points.empty()) {
            Point_2 point = Point_2(q_points.front().x(), q_chain.back().y());
            q_chain.push_back(point);
            q_chain.push_back(q_points.front());
            q_points.pop();
            if (!q_points.empty() && q_chain.back().y() == query_point.y()) {return;}
        }
    }

    void process_segment(std::queue<Point_2>& q_points, const Point_2& source, const Point_2& target, const Point_2& query_point, bool source_in_q, bool target_in_q, bool is_vertical_segment, bool keep_min_x, bool keep_min_y) {
        if (source_in_q && target_in_q && is_vertical_segment) {
            q_points.push(Point_2(source.x(), keep_min_y ? std::min(source.y(), target.y()) : std::max(source.y(), target.y())));
        } else if (source_in_q && target_in_q) {
            q_points.push(Point_2(keep_min_x ? std::min(source.x(), target.x()) : std::max(source.x(), target.x()), source.y()));
        } else if ((source_in_q || target_in_q) && is_vertical_segment) {
            q_points.push(Point_2(source.x(), query_point.y()));
        } else if (source_in_q || target_in_q) {
            q_points.push(Point_2(query_point.x(), source.y()));
        }
    }

    void sort_and_filter_points(std::queue<Point_2>& points_queue, bool ascending_x, bool ascending_y) {
        if (points_queue.empty()) {
            return;
        }

        std::vector<Point_2> points;
        while (!points_queue.empty()) {
            points.push_back(points_queue.front());
            points_queue.pop();
        }

        // Sort by y-coord; if equal by x-coord
        std::sort(points.begin(), points.end(), [ascending_x, ascending_y](const Point_2& a, const Point_2& b) {
            if (a.x() == b.x()) {
                return ascending_y ? a.y() > b.y() : a.y() < b.y();
            }
            return ascending_x ? a.x() < b.x() : a.x() > b.x();
        });

        std::queue<Point_2> all_points;

        for (const auto& point : points) {
            all_points.push(point);
        }

        Point_2 first_point = all_points.front();
        auto current_x = first_point.x();
        auto current_y = first_point.y();
        points_queue.push(first_point);
        all_points.pop();

        while (!all_points.empty()) {
            Point_2 p = all_points.front();
            if (p.x() != current_x && (ascending_y ? p.y() > current_y : p.y() < current_y)) {
                points_queue.push(p);
                current_x = p.x();
                current_y = p.y();
            }
            all_points.pop();
        }
    }

    bool is_vertical(const Segment_2& edge) {
        return edge.source().x() == edge.target().x();
    }
};

class VisibilityPolygonCalculator {
public:
  VisibilityPolygonCalculator(Polygon2WithHoles &poly) {
    if (!poly.outer_boundary().is_simple()) {
      throw std::runtime_error("Polygon is not simple");
    }
    polygon = poly;
    std::vector<Segment2> segments;
    if (poly.outer_boundary().area() < 0) {
      throw std::runtime_error("Polygon is not counterclockwise oriented");
    }
    for (const auto e : poly.outer_boundary().edges()) {
      auto s = e.source();
      auto t = e.target();
      if (s.x() != t.x() && s.y() != t.y()) {
          throw std::runtime_error("Polygon is not rectilinear");
      }
      auto seg = Segment2(s, t);
      segments.push_back(seg);
    }
    for (const auto &hole : poly.holes()) {
      if (hole.area() > 0) {
        throw std::runtime_error("Hole is not clockwise oriented");
      }
      for (const auto e : hole.edges()) {
        auto s = e.source();
        auto t = e.target();
          if (s.x() != t.x() && s.y() != t.y()) {
              throw std::runtime_error("Polygon is not rectilinear");
          }
        auto seg = Segment2(s, t);
        segments.push_back(seg);
      }
    }
    CGAL::insert_non_intersecting_curves(env, segments.begin(), segments.end());
    pl = PointLocation(env);
    auto face = env.unbounded_face();
    if (face->number_of_holes() != 1 || !face->is_unbounded()) {
      throw std::runtime_error("Bad arrangement. Could not determine polygon face.");
    }
    auto hole_it = face->holes_begin();
    assert(hole_it != face->holes_end());
    auto f = (*hole_it)->twin()->face();
    if (f->is_unbounded()) {
      throw std::runtime_error("Bad arrangement. Face should not be unbounded.");
    }
    interior_face = f;
  }

  bool is_feasible_query_point(const Point &query_point) {
    if (polygon.outer_boundary().bounded_side(query_point) ==
        CGAL::ON_UNBOUNDED_SIDE) {
      return false;
    }
    for (const auto &hole : polygon.holes()) {
      if (hole.bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE) {
        return false;
      }
    }
    return true;
  }

  Polygon2 get_visibility_face_handle(const Point &query_point) {
    auto location = pl.locate(query_point);
    const Arrangement_2::Vertex_const_handle *v;
    Rectilinear_Visibility rv(env);
    if ((v = std::get_if<Arrangement_2::Vertex_const_handle>(&location))) {
      auto e_ = (*v)->incident_halfedges();
      if (e_->face() == interior_face) {
          Halfedge_const_handle hh = e_;
        return rv.compute_visibility(query_point, hh);
      }
      ++e_;
      while (e_ != (*v)->incident_halfedges()) {
        if (e_->face() == interior_face) {
            Halfedge_const_handle hh = e_;
          return rv.compute_visibility(query_point, hh);
        }
        ++e_;
      }
      throw std::runtime_error(
          "Query point is on a vertex, but not in interior face");
    }
  }

  Polygon2 compute_visibility_polygon(const Point &query_point) {
    if (!is_feasible_query_point(query_point)) {
      throw std::runtime_error("Query point not feasible");
    }
    return get_visibility_face_handle(query_point);
  }

  Polygon2WithHoles polygon;
  Arrangement_2 env;
  Face_handle interior_face;
  PointLocation pl;
};

template <typename EdgeHandle>
Polygon2 _boundary_to_polygon(const EdgeHandle &e) {
  Polygon2 poly;
  std::vector<Point> points;
  auto e_ = e;
  points.push_back(e_->source()->point());
  ++e_;
  while (e_ != e) {
    points.push_back(e_->source()->point());
    ++e_;
  }
  return Polygon2(points.begin(), points.end());
}

template <typename FaceHandle> Polygon2 _face_to_polygon(const FaceHandle &fh) {
  assert(!fh->is_unbounded());
  return _boundary_to_polygon(fh->outer_ccb());
}

Arrangement_2 _polygon_to_arrangement(const Polygon2WithHoles &poly) {
  // Create a new arrangement from a polygon with holes
  Arrangement_2 env;
  std::vector<Segment2> segments;
  // Add the outer boundary
  for (const auto e : poly.outer_boundary().edges()) {
    auto s = e.source();
    auto t = e.target();
    auto seg = Segment2(s, t);
    segments.push_back(seg);
  }
  // Add the holes
  for (const auto &hole : poly.holes()) {
    for (const auto e : hole.edges()) {
      auto s = e.source();
      auto t = e.target();
      auto seg = Segment2(s, t);
      segments.push_back(seg);
    }
  }
  // Insert the segments into the arrangement
  CGAL::insert(env, segments.begin(), segments.end());
  return env;
}

Ex_arrangement _polygon_to_ex_arrangement(const Polygon2WithHoles &poly) {
  // similar to above; but for Ex_arrangement
  Ex_arrangement env;
  std::vector<Segment2> segments;
  for (const auto e : poly.outer_boundary().edges()) {
    auto s = e.source();
    auto t = e.target();
    auto seg = Segment2(s, t);
    segments.push_back(seg);
  }
  for (const auto &hole : poly.holes()) {
    for (const auto e : hole.edges()) {
      auto s = e.source();
      auto t = e.target();
      auto seg = Segment2(s, t);
      segments.push_back(seg);
    }
  }
  CGAL::insert(env, segments.begin(), segments.end());
  return env;
}

std::vector<Polygon2WithHoles> repair(const Polygon2WithHoles &poly) {
  // Repair a polygon with holes that is self intersecting.
  // Use arrangements to separate the polygons.

  // Create an arrangement
  Arrangement_2 env = _polygon_to_arrangement(poly);
  std::vector<Polygon2WithHoles> result;
  // Get the faces
  for (auto f = env.faces_begin(); f != env.faces_end(); ++f) {
    // Face is a polygon if it is adjacent to the unbounded face
    if (f->is_unbounded()) {
      continue;
    }
    if (!f->outer_ccb()->twin()->face()->is_unbounded()) {
      continue;
    }
    // face to polygon with holes
    // outer boundary
    auto outer_boundary = _face_to_polygon(f);
    // holes
    std::vector<Polygon2> holes;
    for (auto h = f->holes_begin(); h != f->holes_end(); ++h) {
      // h is Ccb_halfedge_circulator
      auto hole_poly = _boundary_to_polygon(*h);
      hole_poly.reverse_orientation();
      holes.push_back(hole_poly);
      assert(holes.back().area() < 0);
    }
    result.push_back(
        Polygon2WithHoles(outer_boundary, holes.begin(), holes.end()));
  }
  return result;
}

struct Guard_overlay {
    std::set<int> operator()(const std::set<int>& set1, const std::set<int>& set2) const {
        std::set<int> result;
        std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(result, result.begin()));
        return result;
    }
};

using Face_overlay_traits = CGAL::Arr_face_overlay_traits<Ex_arrangement, Ex_arrangement, Ex_arrangement, Guard_overlay>;

PYBIND11_MODULE(_cgal_bindings, m) {
  namespace py = pybind11;
  m.doc() = "r-visibility in Python.";

  // Exact numbers
  py::class_<Kernel::FT>(m, "FieldNumber",
                         "A container for exact numbers in CGAL.")
      .def(py::init<long>())
      .def(py::init<double>())
      .def(py::self / Kernel::FT())
      .def(py::self + Kernel::FT())
      .def(py::self - Kernel::FT())
      .def(py::self * Kernel::FT())
      .def(py::self == Kernel::FT())
      .def(py::self < Kernel::FT())
      .def(py::self > Kernel::FT())
      .def(py::self <= Kernel::FT())
      .def(py::self >= Kernel::FT())
      .def("__float__", &CGAL::to_double<Kernel::FT>)
      .def("__str__", [](const Kernel::FT &x) {
        return std::to_string(CGAL::to_double(x));
      });

  // Points
  py::class_<Point>(m, "Point", "A 2-dimensional point")
      .def(py::init<long, long>())
      .def(py::init<double, double>())
      .def(py::init<Kernel::FT, Kernel::FT>())
      .def("x", [](const Point &p) { return p.x(); })
      .def("y", [](const Point &p) { return p.y(); })
      .def("__len__", [](const Point &self) { return 2; })
      .def("__item__",
           [](const Point &self, int i) {
             if (i == 0) {
               return self.x();
             } else if (i == 1) {
               return self.y();
             }
             throw std::out_of_range("Only 0=x and 1=y.");
           })
      .def(py::self == Point())
      .def("__str__", [](const Point &p) {
        return fmt::format("({}, {})", CGAL::to_double(p.x()),
                           CGAL::to_double(p.y()));
      });

  // Polygons
  py::class_<Polygon2>(m, "Polygon", "A simple polygon in CGAL.")
      .def(py::init<>())
      .def(py::init([](const std::vector<Point> &vertices) {
        return std::make_unique<Polygon2>(vertices.begin(), vertices.end());
      }))
      .def("boundary",
           [](const Polygon2 &poly) {
             std::vector<Point> points;
             std::copy(poly.begin(), poly.end(), std::back_inserter(points));
             return points;
           })
      .def("is_simple", &Polygon2::is_simple)
      .def("contains",
           [](const Polygon2 &self, const Point &p) {
             return self.bounded_side(p) != CGAL::ON_UNBOUNDED_SIDE;
           })
      .def("on_boundary",
           [](const Polygon2 &self, const Point &p) {
             return self.bounded_side(p) == CGAL::ON_BOUNDARY;
           })
      .def("area", [](const Polygon2 &poly) { return poly.area(); });

  py::class_<Polygon2WithHoles>(m, "PolygonWithHoles",
                                "A polygon with holes in CGAL.")
      .def(py::init([](const Polygon2 &outer,
                       const std::vector<Polygon2> &holes) {
        for (const auto &hole_poly : holes) {
          if (hole_poly.area() >= 0) {
            throw std::runtime_error("Hole is not clockwise oriented");
          }
        }
        if (outer.area() <= 0) {
          throw std::runtime_error("Polygon is not counterclockwise oriented");
        }
        return Polygon2WithHoles(outer, holes.begin(), holes.end());
      }))
      .def(py::init([](const Polygon2 &outer) {
        if (outer.area() <= 0) {
          throw std::runtime_error("Polygon is not counterclockwise oriented");
        }
        return Polygon2WithHoles(outer);
      }))
      .def(py::init([](const std::vector<Point> &outer_vertices) {
        auto poly = Polygon2(outer_vertices.begin(), outer_vertices.end());
        return Polygon2WithHoles(poly);
      }))
      .def(py::init([](const std::vector<Point> &outer_vertices,
                       const std::vector<std::vector<Point>> &hole_vertices) {
        auto poly = Polygon2(outer_vertices.begin(), outer_vertices.end());
        if (poly.area() <= 0) {
          throw std::runtime_error("Polygon is not counterclockwise oriented");
        }
        std::vector<Polygon2> holes;
        for (auto hole_boundary : hole_vertices) {
          auto hole_poly = Polygon2(hole_boundary.begin(), hole_boundary.end());
          if (hole_poly.area() >= 0) {
            throw std::runtime_error("Hole is not clockwise oriented");
          }
          holes.push_back(hole_poly);
        }
        return Polygon2WithHoles(poly, holes.begin(), holes.end());
      }))
      .def(py::self == Polygon2WithHoles())
      .def(
          "outer_boundary",
          [](const Polygon2WithHoles &self) { return self.outer_boundary(); },
          "Returns a list with all holes (simple polygons).")
      .def("holes",
           [](const Polygon2WithHoles &poly) {
             // Copy the holes into a vector so that PyBind11 can return it
             // as a Python list.
             std::vector<Polygon2> holes;
             std::copy(poly.holes_begin(), poly.holes_end(),
                       std::back_inserter(holes));
             return holes;
           })
      .def("contains",
           [](Polygon2WithHoles &self, const Point &p) {
             if (self.outer_boundary().bounded_side(p) ==
                 CGAL::ON_UNBOUNDED_SIDE) {
               return false;
             }
             for (auto hole : self.holes()) {
               if (hole.bounded_side(p) == CGAL::ON_BOUNDED_SIDE) {
                 return false;
               }
             }
             return true;
           })
      .def("on_boundary",
           [](Polygon2WithHoles &self, const Point &p) {
             if (self.outer_boundary().bounded_side(p) == CGAL::ON_BOUNDARY) {
               return true;
             }
             for (const auto &hole : self.holes()) {
               if (hole.bounded_side(p) == CGAL::ON_BOUNDARY) {
                 return true;
               }
             }
             return false;
           })
      .def("area",
           [](const Polygon2WithHoles &poly) {
             auto area = poly.outer_boundary().area();
             for (const auto &hole : poly.holes()) {
               area -= hole.area();
             }
             return CGAL::to_double(area);
           })
      .def(
          "join",
          [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
            std::vector<Polygon2WithHoles> result;
            Polygon2WithHoles joined;
            if (CGAL::join(self, other, joined,
                           /*UsePolylines=*/CGAL::Tag_false{})) {
              result.push_back(joined);
            } else {
              result.push_back(self);
              result.push_back(other);
            }
            return result;
          },
          "Joins two polygons (with holes) into a list of polygons (with "
          "holes).")
      .def(
          "intersection",
          [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
            std::vector<Polygon2WithHoles> result;
            CGAL::intersection(self, other, std::back_inserter(result),
                               /*UsePolylines=*/CGAL::Tag_false{});
            return result;
          },
          "Computes the intersection of two polygons (with holes). Returns a "
          "list of polygons (with holes).")
      .def(
        "do_intersect",
        [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
          return CGAL::do_intersect(self, other);
        },
        "Check if two polygons (with holes) intersect. Returns a boolean.")
      .def(
          "difference",
          [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
            std::vector<Polygon2WithHoles> result;
            CGAL::difference(self, other, std::back_inserter(result),
                             /*UsePolylines=*/CGAL::Tag_false{});
            return result;
          },
          "Removes the area of the other polygon. Returns a list of polygons "
          "(with holes).");

  py::class_<VisibilityPolygonCalculator>(
      m, "VisibilityPolygonCalculator",
      "A class to compute visibility polygons.")
      .def(py::init<Polygon2WithHoles &>())
      .def("compute_visibility_polygon",
           &VisibilityPolygonCalculator::compute_visibility_polygon,
           py::arg("query_point"),
           "Compute the visibility polygon for a query point.")
      .def("is_feasible_query_point",
           &VisibilityPolygonCalculator::is_feasible_query_point,
           py::arg("query_point"),
           "Check if the query point is within the polygon.");

  py::class_<Ex_arrangement>(m, "AVP_Arrangement",
                            "A class to represent a guard arrangement.")
      .def(py::init<>([](const Polygon2WithHoles &poly, std::set<int> guard) {
            Ex_arrangement arr = _polygon_to_ex_arrangement(poly);
            for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
              if (f->is_unbounded())
                continue;
              f->set_data(guard);
            }
            return arr;
          }))
      .def("overlay",
           [](Ex_arrangement &self, const Ex_arrangement &other) {
             Ex_arrangement result;
             Face_overlay_traits traits;
             CGAL::overlay(self, other, result, traits);
             return result;
           },
           "Computes the overlay of two arrangements.")

      .def("get_shadow_witnesses",
          [](Ex_arrangement &self) {
            std::map<int, std::set<int>> witness_to_guards;
            int index = 1;

            for (auto f = self.faces_begin(); f != self.faces_end(); ++f) {
              if (f->data().empty())
                continue;
              bool is_shadow = true;
              for (auto half_edge = f->outer_ccbs_begin(); half_edge != f->outer_ccbs_end(); ++half_edge) {
                Ex_arrangement::Ccb_halfedge_circulator curr = *half_edge;
                Ex_arrangement::Ccb_halfedge_circulator done = curr;
                do {
                  if (!is_shadow)
                    break;
                  if (curr->twin()->face()->data().empty()) {
                    ++curr;
                    continue;
                  }
                  if (std::includes(f->data().begin(), f->data().end(), curr->twin()->face()->data().begin(), curr->twin()->face()->data().end())) {
                    is_shadow = false;
                  }
                  ++curr;
                } while (curr != done);
              }

              if (is_shadow) {
                witness_to_guards.insert({index, f->data()});
                index++;
              }
            }
            return witness_to_guards;
          },
          "Returns a dictionary with shadow witness and corresponding guard positions.");

  py::class_<Arrangement_2>(m, "Arrangement",
                            "Class for CGAL::Arrangement_2.")
      .def(py::init<>([](const Polygon2WithHoles &poly) {
            return _polygon_to_arrangement(poly);
          }))
      .def("overlay",
           [](Arrangement_2 &self, const Arrangement_2 &other) {
             Arrangement_2 result;
             CGAL::overlay(self, other, result);
             return result;
           },
           "Computes the overlay of two arrangements.");

  m.def("repair", &repair,
        "Repair a polygon with holes that is self "
        "intersecting. Returns a list of polygons with "
        "holes.");
}
