#include "Generator.h"
#include <iostream>

//using namespace CGAL::parameters;
//
//// Domain
//typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
//typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
//        Function_wrapper;
//typedef Function_wrapper::Function_vector Function_vector;
//typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
//// Triangulation
//typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
//typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
//// Mesh Criteria
//typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
//typedef Mesh_criteria::Facet_criteria    Facet_criteria;
//typedef Mesh_criteria::Cell_criteria     Cell_criteria;
//
//typedef K::FT FT;
//typedef K::Point_3 Point;
//
//Point center_1(0,0,0);

//std::function<int (int)> retFun() {
//    return [](int x) { return x; };
//}

int fonc() {
//
//    // Define functions
//    Function f1(para(0.5,0.5,0.5));
//    Function f2(&sphere_function);
//    Function f3(&torus_function);
//    Function f4(&cube_function);
//    Function_vector v;
//    v.push_back(f1);
//    v.push_back(f3);
//    std::vector<std::string> vps;
//    vps.push_back("-+");
//
//    namespace param = CGAL::parameters;
//    Mesh_domain domain(param::function = Function_wrapper(v, vps),
//                       param::bounding_object = K::Sphere_3(CGAL::ORIGIN,
//                                                            5.*5.),
//                       param::relative_error_bound = 1e-6);
//    // Set mesh criteria
//    Facet_criteria facet_criteria(30, 0.2, 0.02); // angle, size, approximation
//    Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
//    Mesh_criteria criteria(facet_criteria, cell_criteria);
//
//    // Mesh generation
//    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::no_exude(), CGAL::parameters::no_perturb());
//    // Output
//    std::ofstream medit_file("mymesh");
//    c3t3.output_to_medit(medit_file);
//    // Mesh_domain domain =
//    //         Mesh_domain::create_implicit_mesh_domain(sphere_function,
//    //                                                  K::Sphere_3(CGAL::ORIGIN, 2.));
//    // // Mesh criteria
//    // Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
//    //                        cell_radius_edge_ratio=2, cell_size=0.4);
//    //
//    //  // Mesh generation
//    // C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
//    //  // Output
//    // for (auto it =  c3t3.; it != c3t3.vertices_in_complex_end(); it++){
//    //     c3t3.
//    //     std::cout << "ola";
//    // }
//
//    int tot = 0;
//
//    Tr& tr = c3t3.triangulation();
//
//    for( Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
//         vit != tr.finite_vertices_end();
//         ++vit)
//    {
//        Tr::Weighted_point p = tr.point(vit);
//        std::cout << CGAL::squared_distance(Point(center_1), p) << '\n';
//    //           << get(vertex_pmap, vit)
//    //           << '\n';
//        tot++;
//    }
//
//    //    std::cout << c3t3.REGULAR
//
//    //    c3t3.cor
//    //    c3t3.output_to_medit(medit_file);

    return 0;
}

using std::vector;

void test(){
    vector<double> dimens;
    dimens.push_back(1.0);
    dimens.push_back(1.0);
    dimens.push_back(1.0);

    vector<double> c1 = {0.1, 0.5};
    vector<double> c2 = {0.1, 0.5};
    vector<double> r = {0.2, 0.1};
    vector<int> a = {1, 1};
    

    vector<MeshNode>nodes = generateMesh(dimens, c1, c2, r, a, 1, 1, 1, 1, 1);

    using std::cout;
    for (MeshNode node : nodes){
        cout << "Node index: " << node.index << "\n";
        cout << "x: " << node.x << " y: " << node.y << " z: " << node.z << "\n";
        cout << "Adjacency: ";
        for (auto ind: node.adj)
            cout << ind << " ";
        cout << "\n";
    }
}

int main(){
    test();
}
