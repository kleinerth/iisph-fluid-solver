#include "io/parameterIO.h"
#include "fluid/particleBasics.h"
#include <fstream> // std::ifstream, std::ofstream
#include <iostream>
#include <sstream> // std::stringstream
#include <string>

namespace iisph::io {

SimulationParameters loadSimulationParameters(const std::string& filename) {
  SimulationParameters params{};

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open parameter file: " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string key, value;

    if (std::getline(iss, key, '=') && std::getline(iss, value)) {
      std::istringstream valStream(value);

      if (key == "startFromSavedState")
        valStream >> params.startFromSavedState;
      else if (key == "numIterations")
        valStream >> params.numIterations;
      else if (key == "snapshotFrequency")
        valStream >> params.snapshotFrequency;
      else if (key == "saveState")
        valStream >> params.saveState;
      else if (key == "maxJacobiIterations")
        valStream >> params.maxJacobiIterations;
      else if (key == "h")
        valStream >> params.h;
      else if (key == "restDensity")
        valStream >> params.restDensity;
      else if (key == "ny")
        valStream >> params.ny;
      else if (key == "deltaT")
        valStream >> params.deltaT;
      else if (key == "toleratedError")
        valStream >> params.toleratedError;
      else if (key == "gamma")
        valStream >> params.gamma;
      else if (key == "gammaDens")
        valStream >> params.gammaDens;
      else if (key == "zeroTolerance")
        valStream >> params.zeroTolerance;
      else if (key == "gravity")
        valStream >> params.gravity;
      else if (key == "omega")
        valStream >> params.omega;
      else if (key == "CFL")
        valStream >> params.CFL;
      else if (key == "baseLength")
        valStream >> params.baseLength;
      else if (key == "wallHeight")
        valStream >> params.wallHeight;
      else if (key == "wallThickness")
        valStream >> params.wallThickness;
      else if (key == "fluidHeight")
        valStream >> params.fluidHeight;
      else if (key == "fluidWidth")
        valStream >> params.fluidWidth;
    }
  }
  return params;
}

} // namespace iisph::io
