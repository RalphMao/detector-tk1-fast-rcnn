#include <boost/python.hpp>
#include <numpy/arrayobject.h> 
#include <Python.h> 

namespace bp = boost::python;

bp::object feed_out(bp::object input)
{
    PyArrayObject* array= reinterpret_cast<PyArrayObject *>(input.ptr());

    bp::object* return_array = reinterpret_cast<bp::object *>(array);
    return *return_array;
}

char const* greet()
{
       return "hello, world";
}
BOOST_PYTHON_MODULE(hello_ext)
{
        using namespace boost::python;
            def("greet", greet);
            def("feed_out",feed_out);
}
