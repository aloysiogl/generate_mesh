#include<iostream>
#include<vector>
#include <CGAL/Point_2.h>
#include <CGAL/Cartesian.h>
#include "Test.h"
#include "TestInc.h"

using namespace std;

void print_something(){
  CGAL::Point_2< CGAL::Cartesian<double> > p(1.0, 2.0);
  std::cout << p.x();
}
vector<myobj> print_vector(std::vector<int> v){
//    for (int i = 0; i < v.size(); ++i)
//        cout << v[i] << " ";
    myobj a(1);
    vector<myobj> vec;
    vec.push_back(a);
    return vec;
}
