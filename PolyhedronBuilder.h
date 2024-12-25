//#include <vector>
//#include <fstream>
//#include "struct.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
//#include <CGAL/Mesh_criteria_3.h>
//#include <CGAL/make_mesh_3.h>
//#include <CGAL/Mesh_polyhedron_3.h>
//#include <CGAL/Polyhedron_incremental_builder_3.h>
//#include <CGAL/Mesh_domain_with_polyline_features_3.h>
//#include <CGAL/IO/Polyhedron_iostream.h>
//
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
//typedef CGAL::Mesh_polyhedron_3<Polyhedron>::type Mesh_polyhedron;
//typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> Mesh_domain3;
//typedef Mesh_domain3::Index Index;
//typedef Mesh_domain3::Corner_index Corner_index;
//typedef Mesh_domain3::Surface_patch_index Facet_index;
//typedef Mesh_domain3::Index Halfedge_index;
//typedef Mesh_domain3::Index Vertex_index;
//typedef Kernel::Point_3 Point_3;
//typedef Polyhedron::HalfedgeDS HalfedgeDS;
//
//typedef CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain> Mesh_domain_with_features;
//
//class PolyhedronBuilder : public CGAL::Modifier_base<HalfedgeDS> {
//public:
//    PolyhedronBuilder(std::vector<Point>& points) : m_points(points), m_triangle_index(0) {}
//
//    void operator()(HalfedgeDS& hds) {
//        typedef typename HalfedgeDS::Vertex   Vertex;
//        typedef typename Vertex::Point       Point;
//        CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(hds, true);
//        builder.begin_surface(m_points.size(), m_points.size() * 3);
//
//        for (const auto& point : m_points) {
//            const Point p(point.x, point.y, point.z);
//            builder.add_vertex(p);
//        }
//
//        for (std::size_t i = 0; i < m_points.size(); i += 3) {
//            builder.begin_facet();
//            builder.add_vertex_to_facet(i);
//            builder.add_vertex_to_facet(i + 1);
//            builder.add_vertex_to_facet(i + 2);
//            builder.end_facet();
//        }
//
//        builder.end_surface();
//    }
//
//    Polyhedron build_polyhedron(std::vector<Point>& points) {
//        Polyhedron P;
//        PolyhedronBuilder builder(points);
//        P.delegate(builder);
//        return P;
//    }
//
//
//    void write_polyhedron_obj(const Polyhedron& polyhedron, const std::string& filename) {
//        std::ofstream out_file(filename);
//    }
//
//private:
//    std::vector<Point>& m_points; 
//    int m_triangle_index;
//};
