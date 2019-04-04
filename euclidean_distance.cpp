#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <math.h>

double metric(double x1, double y1, double x2, double y2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

PYBIND11_MODULE(euclidean_distance, m) {
    m.doc() = "euclidean distance metric";
    m.def("metric", &metric, "A function which returns the euclidean distance from (x1,y1) to (x2,y2)");
}
