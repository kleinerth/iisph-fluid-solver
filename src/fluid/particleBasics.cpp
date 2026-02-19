#include "fluid/particleBasics.h"
#include <cmath>  // sqrt
#include <limits> // numeric limits

namespace iisph::fluid {

BoundingBox computeBoundingBox(const Particle* particles, size_t numParticles) {
  double xmin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < numParticles; i++) {
    const Particle& p = particles[i];
    if (p.posX < xmin) {
      xmin = p.posX;
    }
    if (p.posY < ymin) {
      ymin = p.posY;
    }
    if (p.posX > xmax) {
      xmax = p.posX;
    }
    if (p.posY > ymax) {
      ymax = p.posY;
    }
  }

  return BoundingBox{xmin, xmax, ymin, ymax};
}

double euclideanDistance(double x1, double y1, double x2, double y2) {
  return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

} // namespace iisph::fluid
