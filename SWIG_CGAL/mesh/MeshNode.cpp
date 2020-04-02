//
// Created by Aloysio Galvão Lopes on 2020-03-31.
//

#include "MeshNode.h"

MeshNode::MeshNode() : x(0.0), y(0.0), z(0.0), index(-1){}

MeshNode::MeshNode(int index) : x(0.0), y(0.0), z(0.0), index(index) {}

MeshNode::MeshNode(double x, double y, double z, int index) : x(x), y(y), z(z), index(index) {}

MeshNode::~MeshNode() {
    
}

