#ifndef PARAMETERIO_H
#define PARAMETERIO_H

#include "fluid/particleBasics.h"
#include <string>

namespace iisph::io {

struct SimulationParameters {
  bool startFromSavedState;

  int numIterations;
  int snapshotFrequency;
  int saveState;
  int maxJacobiIterations;

  double h;
  double restDensity;
  double ny;
  double deltaT;
  double toleratedError;
  double gamma;
  double gammaDens;
  double zeroTolerance;
  double gravity;
  double omega;
  double CFL;

  size_t baseLength;
  size_t wallHeight;
  size_t wallThickness;
  size_t fluidHeight;
  size_t fluidWidth;
};

SimulationParameters loadSimulationParameters(const std::string& filename);

} // namespace iisph::io

#endif
