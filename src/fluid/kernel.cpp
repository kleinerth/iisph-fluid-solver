#include "fluid/kernel.h"
#include "fluid/particleBasics.h"
#include <algorithm> // std::max
#include <utility>   // std::pair

double W(double x1, double y1, double x2, double y2, double h, double alpha) {
  double dist = iisph::fluid::euclideanDistance(x1, y1, x2, y2) / h;
  double t1 = std::max(1.0 - dist, 0.0);
  double t2 = std::max(2.0 - dist, 0.0);
  double kernelValue = alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);
  return kernelValue;
}

std::pair<double, double> gradW(double x1, double y1, double x2, double y2, double h, double alpha) {
  double dist = iisph::fluid::euclideanDistance(x1, y1, x2, y2) / h;
  // Gradient between a particle and itself is zero
  if (dist < 1e-12) {
    return std::make_pair(0.0, 0.0);
  }
  double t1 = std::max(1.0 - dist, 0.0);
  double t2 = std::max(2.0 - dist, 0.0);
  std::pair<double, double> gradient;
  gradient.first = alpha * ((x1 - x2) / (dist * h * h)) * ((-3) * t2 * t2 + 12 * t1 * t1);
  gradient.second = alpha * ((y1 - y2) / (dist * h * h)) * ((-3) * t2 * t2 + 12 * t1 * t1);
  return gradient;
}
