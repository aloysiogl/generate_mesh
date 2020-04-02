//
// Created by Aloysio Galvão Lopes on 2020-03-31.
//

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#ifndef MESH_FTTOPOINTFUNCTIONWRAPPER_H
#define MESH_FTTOPOINTFUNCTIONWRAPPER_H

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

template <typename FT, typename P>
class FT_to_point_function_wrapper : public CGAL::cpp98::unary_function<P, FT>
{
    typedef FT (*Implicit_function)(FT, FT, FT);
    Implicit_function functionFixedPointer;
    std::function<double(double, double, double)> functionVariablePointer;
    bool useFixedPointer = true;

public:
    typedef P Point;

    FT_to_point_function_wrapper(Implicit_function function): functionFixedPointer(function) {
        // Just initialize the function
    }
    FT_to_point_function_wrapper(std::function<double(double, double, double)> f) {
        // Initialize the function class to be able to handle run time
        // generated functions
        functionVariablePointer = f;
        useFixedPointer = false;
    }

    FT operator()(Point p) const {
        // Selecting the right implementation of the function
        if (useFixedPointer) {
            return functionFixedPointer(p.x(), p.y(), p.z());
        }
        return functionVariablePointer(p.x(), p.y(), p.z());
    }
};


#endif //MESH_FTTOPOINTFUNCTIONWRAPPER_H
