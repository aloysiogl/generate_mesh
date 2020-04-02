//
// Created by Aloysio Galvão Lopes on 2020-03-31.
//

#include <functional>

#ifndef MESH_IMPLICITFUNCTIONS_H
#define MESH_IMPLICITFUNCTIONS_H

typedef std::function<double(double, double, double)> fptr;

fptr parallelepiped_function (double xMax, double yMax, double zMax){
    auto lamb = std::function<double(double, double, double)>([=](double x, double y, double z) {
        if( x < xMax && x > 0 &&
            y < yMax && y > 0 &&
            z < zMax && z > 0 )
            return -1.;
        return 1.;
    });
    return lamb;
}

fptr cylinder_function (double c1, double c2, double r, int axis) {
    auto lamb = std::function<double(double, double, double)>([=](double x, double y, double z) {
        if (axis == 0){
            double y2=(y-c1)*(y-c1), z2=(z-c2)*(z-c2);
            return y2+z2-r*r;
        } else if (1 == axis){
            double x2=(x-c1)*(x-c1), z2=(z-c2)*(z-c2);
            return x2+z2-r*r;
        } else {
            double x2=(x-c1)*(x-c1), y2=(y-c2)*(y-c2);
            return x2+y2-r*r;
        }
    });
    return lamb;
}

#endif //MESH_IMPLICITFUNCTIONS_H
