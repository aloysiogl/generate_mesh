%module (package="CGAL") mesh

//include files
%{
    #include <SWIG_CGAL/mesh/Generator.h>
    #include <SWIG_CGAL/mesh/MeshNode.h>
    #include <SWIG_CGAL/mesh/mytest.h>
    #include <SWIG_CGAL/mesh/Test.h>
    #include <SWIG_CGAL/mesh/TestInc.h>
%}


%include "std_vector.i"
// Instantiate templates used by example
namespace std {
        %template(IntVector) vector<int>;
        %template(DoubleVector) vector<double>;
        %template(myobjVector) vector<myobj>;
        %template(MeshNodeVector) vector<MeshNode>;
}

%include "SWIG_CGAL/mesh/Generator.h"
%include "SWIG_CGAL/mesh/MeshNode.h"
%include "SWIG_CGAL/mesh/mytest.h"
%include "SWIG_CGAL/mesh/Test.h"
%include "SWIG_CGAL/mesh/TestInc.h"