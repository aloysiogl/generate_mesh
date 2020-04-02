//
// Created by Aloysio Galvão Lopes on 2020-03-31.
//

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include "implicit_functions.h"
#include "FTtoPointFunctionWrapper.h"

#include "FTtoPointFunctionWrapper.h"
#include "ImplicitFunctions.h"
#include "Generator.h"

// Domain
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
        Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

typedef K::FT FT;
typedef K::Point_3 Point;

namespace param = CGAL::parameters;

struct point {
    double x, y, z;
    point() {}
    point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    double norm() { return hypot(x, y); }

    point operator +(point p) const {
        return point(x + p.x, y + p.y, z + p.z);
    }
    point operator -(point p) const {
        return point(x - p.x, y - p.y, z - p.z);
    }
    point operator *(double k) const {
        return point(k*x, k*y, k*z);
    }
    point operator /(double k) const {
        return point(x/k, y/k, z/k);
    }
    double dist(point p) const {
        return hypot(hypot(x-p.x, y-p.y), z-p.z);
    }
    bool operator ==(point p) const {
        return dist(p) < 1e-6;
    }
    double inner(point p) {
        return x*p.x + y*p.y + z*p.z;
    }
};

double distanceSegs(vector<double> &dimensions, double c1, double c2, int axis, point &p3, point &p4) {
    // Calculating cylinder segment
    point p1, p2;
    if (axis == 0) {
        p1.x = 0;
        p2.x = dimensions[0];
        p1.y = c1;
        p2.y = c1;
        p1.z = c2;
        p2.z = c2;
    } else if (axis == 1) {
        p1.x = c1;
        p2.x = c1;
        p1.y = 0;
        p2.y = dimensions[1];
        p1.z = c2;
        p2.z = c2;
    } else {
        p1.x = c1;
        p2.x = c1;
        p1.y = c2;
        p2.y = c2;
        p1.z = 0;
        p2.z = dimensions[2];
    }

    point u = p2 - p1;
    point v = p4 - p3;
    point w = p1 - p3;
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
    std::cout << sc << " " << tc << std::endl;

    // get the difference of the two closest points
    point dP = w + (u * sc) - (v * tc);  // =  S1(sc) - S2(tc)
    std::cout << dP.x << " " << dP.y << " " << dP.z << std::endl;

    return std::sqrt(dP.inner(dP));   // return the closest distance
}

vector<MeshNode> generateMesh(vector<double> dimensions, vector<double> positions_1,
                              vector<double> positions_2, vector<double> radius, vector<int> axis) {
    // TODO consider using const references
    assert(positions_1.size() == positions_2.size() && positions_1.size() == radius.size()
           && positions_1.size() == axis.size() && dimensions.size() == 3);

    // Getting function for room generation
    Function room(parallelepiped_function(
            dimensions[0],
            dimensions[1],
            dimensions[2]
            ));

    // Adding functions together
    Function_vector fv;
    vector<string> domainOperations;

    // Adding room mesh
    fv.push_back(room);
    domainOperations.emplace_back("-");

    // Adding cylinders
    for (int i = 0; i < positions_1.size(); ++i){
        Function cyl(cylinder_function(positions_1[i], positions_2[i],
                radius[i], axis[i]));
        fv.push_back(cyl);
        domainOperations[0]+= "+";
    }

    // Creating mesh domain
    // TODO fix the radius
    Mesh_domain domain(param::function = Function_wrapper(fv, domainOperations),
                       param::bounding_object = K::Sphere_3(CGAL::ORIGIN,
                                                            5.*5.),
                       param::relative_error_bound = 1e-6);

    // Set mesh criteria
    // TODO create a file with the criteria
    Facet_criteria facet_criteria(30, 0.2, 0.02); // angle, size, approximation
    Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::no_exude(), CGAL::parameters::no_perturb());

    // Mesh vertices list generation
    vector<MeshNode> nodes;

    // Getting the triangulation with edges and vertices
    Tr& tr = c3t3.triangulation();

    // Adding nodes to output vector
    int nNodes = 0;
    boost::unordered_map<Tr::Vertex_handle, int> verticeHandleToIndice;
    for( Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end();
         ++vit)
    {
        verticeHandleToIndice[vit] = nNodes+1;
        Tr::Weighted_point p = tr.point(vit);
        nodes.emplace_back(p.x(), p.y(), p.z(), nNodes);
        nNodes++;
    }

    // Adding adjacency list to nodes
    for( auto eit = tr.edges_begin() ;
         eit != tr.edges_end() ;
         ++eit )
    {
        // Getting handles to nodes
        Tr::Vertex_handle vh1 = eit->first->vertex((eit->second)%4);
        Tr::Vertex_handle vh2 = eit->first->vertex((eit->third)%4);

        // Getting integer indices
        int i1 = verticeHandleToIndice[vh1];
        int i2 = verticeHandleToIndice[vh2];

        // Verify intersection
        bool validEdge = true;
        for (int i = 0; i < positions_1.size(); ++i){
            point p1 = point(nodes[i1].x, nodes[i1].y, nodes[i1].z);
            point p2 = point(nodes[i2].x, nodes[i2].y, nodes[i2].z);
            if (distanceSegs(dimensions, positions_1[i], positions_2[i], axis[i], p1, p2) < radius[i]){
                validEdge = false;
                break;
            }
        }

        // Excluding index 0 (invalid) and adding to
        // the adjacency list of each node
        if (i1 > 0 && i2 > 0 && validEdge){
            i1--; i2--;
            nodes[i1].adj.push_back(i2);
            nodes[i2].adj.push_back(i1);
        }

    }

    // TODO remove this and the mesh output
    std::cout << "nver: " << tr.number_of_vertices() << "\n";

//    c3t3.cor
    std::ofstream medit_file("mymesh");
    c3t3.output_to_maya(medit_file);

    return nodes;
}
