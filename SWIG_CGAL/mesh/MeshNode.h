//
// Created by Aloysio Galvão Lopes on 2020-03-31.
//

#include <vector>

#ifndef MESH_MESHNODE_H
#define MESH_MESHNODE_H

using std::vector;

struct MeshNode {
    int index;
    double x, y, z;
    vector<int> adj;

    MeshNode();
    MeshNode(int index);
    MeshNode(double x, double y, double z, int index);

    ~MeshNode();
};

#endif //MESH_MESHNODE_H
