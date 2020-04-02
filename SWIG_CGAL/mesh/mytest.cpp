#include "Generator.h"
#include <iostream>
#include "Point.h"
#include <cfloat>

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

double distanceSegs(Point &p1, Point &p2, Point &p3, Point &p4) {
    Point u = p2 - p1;
    Point v = p4 - p3;
    Point w = p1 - p3;
    double a = u.inner(u);         // always >= 0
    double b = u.inner(v);
    double c = v.inner(v);         // always >= 0
    double d = u.inner(w);
    double e = v.inner(w);
    double D = a*c - b*b;        // always >= 0
    double sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
    double tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < 1e-6) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (std::fabs(sN) < 1e-6 ? 0.0 : sN / sD);
    tc = (std::fabs(tN) < 1e-6 ? 0.0 : tN / tD);

    // get the difference of the two closest points
    Point dP = w + (u * sc) - (v * tc);  // =  S1(sc) - S2(tc)

    return std::sqrt(dP.inner(dP));   // return the closest distance
}

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

/*void test(){
    vector<double> dimens;
    dimens.push_back(1.0);
    dimens.push_back(1.0);
    dimens.push_back(1.0);

    vector<double> c1 = {0.1, 0.5};
    vector<double> c2 = {0.1, 0.5};
    vector<double> r = {0.2, 0.1};
    vector<int> a = {1, 1};
    

    vector<MeshNode>nodes = generateMesh(dimens, c1, c2, r, a);

    using std::cout;
    for (MeshNode node : nodes){
        cout << "Node index: " << node.index << "\n";
        cout << "x: " << node.x << " y: " << node.y << " z: " << node.z << "\n";
        cout << "Adjacency: ";
        for (auto ind: node.adj)
            cout << ind << " ";
        cout << "\n";
    }
}*/

double approxDistanceSegs(Point &p1, Point &p2, Point &p3, Point &p4, int ptsPerSeg = 100) {
    Point d1 = (p2-p1) / ptsPerSeg;
    Point d2 = (p4-p3) / ptsPerSeg;
    double minDist = DBL_MAX;
    for (Point a = p1; !(a == p2); a = a + d1) {
        for (Point b = p3; !(b == p4); b = b + d2) {
            minDist = std::min(minDist, a.dist(b));
        }
    }
    return minDist;
}

Point randPt() {
    double x = (rand() % 100) / 50.0 - 2.0;
    double y = (rand() % 100) / 50.0 - 2.0;
    double z = (rand() % 100) / 50.0 - 2.0;
    return Point(x, y, z);
}

int main(){
//    test();
    srand(time(NULL));
    for (int k = 0; k < 1000; k++) {
        Point p1 = randPt();
        Point p2 = randPt();
        Point p3 = randPt();
        Point p4 = randPt();
        
        if (std::abs(distanceSegs(p1, p2, p3, p4) - approxDistanceSegs(p1, p2, p3, p4)) > 0.02) {
            std::cout << "Seg1: (" << p1.x << ", " << p1.y << ", " << p1.z  
                      << ") -> (" << p2.x << ", " << p2.y << ", " << p2.z << ")" << std::endl;
            std::cout << "Seg2: (" << p3.x << ", " << p3.y << ", " << p3.z
                      << ") -> (" << p4.x << ", " << p4.y << ", " << p4.z << ")" << std::endl;
            std::cout << "Dist: " << distanceSegs(p1, p2, p3, p4) << " approx: "
                      << approxDistanceSegs(p1, p2, p3, p4) << std::endl;
        }
    }
}
