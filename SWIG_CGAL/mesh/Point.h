//
// Created by Aloysio Galvão Lopes on 2020-04-02.
//

#ifndef MESH_POINT_H
#define MESH_POINT_H


#include <cmath>

struct Point {
    double x, y, z;
    Point() {}
    Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    double norm() { return hypot(x, y); }

    Point operator +(Point p) const {
        return {x + p.x, y + p.y, z + p.z};
    }
    Point operator -(Point p) const {
        return {x - p.x, y - p.y, z - p.z};
    }
    Point operator *(double k) const {
        return {k*x, k*y, k*z};
    }
    Point operator /(double k) const {
        return {x/k, y/k, z/k};
    }
    double dist(Point p) const {
        return hypot(hypot(x-p.x, y-p.y), z-p.z);
    }
    bool operator ==(Point p) const {
        return dist(p) < 1e-6;
    }
    double inner(Point p) {
        return x*p.x + y*p.y + z*p.z;
    }
};


#endif //MESH_POINT_H
