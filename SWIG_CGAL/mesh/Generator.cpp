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

        //

        // Excluding index 0 (invalid) and adding to
        // the adjacency list of each node
        if (i1 > 0 && i2 > 0){
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