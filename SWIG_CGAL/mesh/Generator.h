//
// Created by Aloysio Galvão Lopes on 2020-03-31.
//

#include <vector>
#include <string>
#include "MeshNode.h"

#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

using std::vector;
using std::string;

vector<MeshNode> generateMesh(vector<double> dimensions, vector<double> positions_1,
                              vector<double> positions_2, vector<double> radius, vector<int> axis,
                              double criteria_angle, double criteria_size, double 	criteria_approximation, double criteria_rediusedge, double criteria_ratio);

#endif //MESH_GENERATOR_H
