#ifndef PARTICLEBASICS_H
#define PARTICLEBASICS_H

#include <cstddef>  // size_t

namespace iisph::fluid {

class Particle {
public:
  double posX, posY, density, pressure, mass, vx, vy, vxPredicted, vyPredicted, ax, ay, sourceTerm, diagElem, densError,
      Ap;
  bool fluid;
  int cell_ID, trackingID;

  Particle()
      : posX(0), posY(0), density(1), pressure(0), mass(1), vx(0), vy(0), vxPredicted(0), vyPredicted(0), ax(0), ay(0),
        sourceTerm(0), diagElem(0), densError(0), Ap(0), fluid(true), cell_ID(0), trackingID(0) {}

  Particle(double posX, double posY, double density, double mass, bool fluid)
      : posX(posX), posY(posY), density(density), pressure(0), mass(mass), vx(0), vy(0), vxPredicted(0), vyPredicted(0),
        ax(0), ay(0), sourceTerm(0), diagElem(0), densError(0), Ap(0), fluid(fluid), cell_ID(0), trackingID(0) {}
};

struct BoundingBox {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

BoundingBox computeBoundingBox(const Particle* particles, size_t numParticles);
double euclideanDistance(double x1, double y1, double x2, double y2);

} // namespace iisph::fluid

#endif
