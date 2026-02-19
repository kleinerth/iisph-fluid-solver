#ifndef KERNEL_H
#define KERNEL_H

#include <utility>   // std::pair

double W(double x1, double y1, double x2, double y2, double h, double alpha);
std::pair<double, double> gradW(double x1, double y1, double x2, double y2, double h, double alpha);

#endif
