cmake_minimum_required(VERSION 3.14)

project(MESH)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE DEBUG)

find_package(CGAL REQUIRED PATHS /Users/aloysiogl/Downloads/CGAL-4.14.3/build)
find_package(Boost REQUIRED COMPONENTS thread)

#add_subdirectory(Source)
#
#add_library(Source SWIG_CGAL/mesh/mytest.cpp)
#
#target_include_directories(
#        Source PUBLIC SWIG_CGAL/mesh
#)

add_executable(mesh SWIG_CGAL/mesh/mytest.cpp
        SWIG_CGAL/mesh/implicit_functions.h
        SWIG_CGAL/mesh/mytest.h
        SWIG_CGAL/mesh/FTtoPointFunctionWrapper.h SWIG_CGAL/mesh/Generator.cpp SWIG_CGAL/mesh/Generator.h SWIG_CGAL/mesh/ImplicitFunctions.h SWIG_CGAL/mesh/MeshNode.cpp SWIG_CGAL/mesh/MeshNode.h SWIG_CGAL/mesh/Point.cpp SWIG_CGAL/mesh/Point.h)
target_link_libraries(mesh  CGAL::CGAL Boost::thread)