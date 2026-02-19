#include "fluid/particleBasics.h"
#include "fluid/kernel.h"
#include "fluid/neighborSearch.h
#include "fluid/setupScene.h"
#include "fluid/simMethods.h"
#include "io/parameterIO.h"
#include "io/particleIO.h"
#include "visualization/graphics.h"
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
  iisph::io::SimulationParameters params = iisph::io::loadSimulationParameters("simParameters.txt");

  // Coefficient for cubic spline kernel in 2D
  const double alpha = 5 / (14 * 3.14159265358979 * params.h * params.h);

  const size_t numParticles = params.wallThickness * params.wallHeight * 2 +
                              params.baseLength * params.wallThickness * 2 + params.fluidHeight * params.fluidWidth;

  iisph::fluid::Particle* particles = new iisph::fluid::Particle[numParticles];

  if (params.startFromSavedState) {
    iisph::io::loadParticlesFromCSV(particles, numParticles, "particle_state");
  } else {
    createDam(particles, numParticles, params.h, params.restDensity, params.baseLength, params.wallHeight,
              params.wallThickness, params.fluidHeight, params.fluidWidth);
  }

  std::vector<int>* neighbors = new std::vector<int>[numParticles];
  for (size_t i = 0; i < numParticles; ++i) {
    neighbors[i].reserve(16);
  }

  for (size_t i = 0; i < numParticles; i++) {
    particles[i].trackingID = i;
    particles[i].vx = 0;
    particles[i].vy = 0;
  }

  sf::RenderWindow window(sf::VideoMode(800, 800), "SPH Simulation");
  float simHeight =
      iisph::visualization::setupWindow(window, params.baseLength, params.h, params.wallHeight, params.wallThickness);

  double avgJacobiIterations = 0.0;
  double errorInPercent = 1.0;
  int itJacobi = 0;

  std::chrono::duration<double> totalNeighborTime(0.0);
  std::chrono::duration<double> totalInitTime(0.0);
  std::chrono::duration<double> totalSolverTime(0.0);

  for (int it = 0; it < params.numIterations; it++) {
    // Compute bounding box around particles for the grid used in neighbor search
    iisph::fluid::BoundingBox box = iisph::fluid::computeBoundingBox(particles, numParticles);

    auto startNeighborSearch = std::chrono::high_resolution_clock::now();
    // First iteration: particles are not sorted, therefore insertion sort has bad average runtime.
    // All iterations except first one: particles are nearly sorted, therefore insertion sort has near linear runtime.
    if (it == 0) {
      iisph::fluid::zIndexSort(particles, neighbors, numParticles, params.h, box.xmin, box.xmax, box.ymin, box.ymax, 2);
    } else {
      iisph::fluid::zIndexSort(particles, neighbors, numParticles, params.h, box.xmin, box.xmax, box.ymin, box.ymax, 1);
    }
    auto endNeighborSearch = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> neighborSearchDuration = endNeighborSearch - startNeighborSearch;
    totalNeighborTime += neighborSearchDuration;

    auto startInit = std::chrono::high_resolution_clock::now();

    iisph::fluid::computeDensity(particles, numParticles, neighbors, params.h, alpha, params.gammaDens);
    iisph::fluid::computeNonPressureAcc(particles, numParticles, neighbors, params.h, alpha, params.ny, params.gravity);
    iisph::fluid::predictVelocity(particles, numParticles, params.deltaT);
    iisph::fluid::computeSourceTerm(particles, numParticles, neighbors, params.h, alpha, params.restDensity,
                                    params.deltaT);
    iisph::fluid::computeDiagonalElement(particles, numParticles, neighbors, params.h, alpha, params.gamma,
                                         params.deltaT);

    for (size_t i = 0; i < numParticles; i++) {
      if (particles[i].fluid) {
        // Save one Jacobi iteration by initializing with the value after the first iteration instead of zero
        particles[i].pressure = std::max(params.omega * particles[i].sourceTerm / particles[i].diagElem, 0.0);
      }
    }

    itJacobi = 0;
    errorInPercent = 1.0;
    double avgDensError = 100.0;

    auto endInit = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> initDuration = endInit - startInit;
    totalInitTime += initDuration;

    auto startSolver = std::chrono::high_resolution_clock::now();

    while (itJacobi < 3 || (itJacobi < params.maxJacobiIterations && errorInPercent * 100 > params.toleratedError)) {
      iisph::fluid::computePressureAcc(particles, numParticles, neighbors, params.h, alpha, params.gamma);
      iisph::fluid::computeAp(particles, numParticles, neighbors, params.h, alpha, params.deltaT);

      for (size_t i = 0; i < numParticles; i++) {
        if (particles[i].fluid) {
          if (std::abs(particles[i].diagElem) > params.zeroTolerance) {
            particles[i].pressure =
                std::max(0.0, particles[i].pressure +
                                  params.omega * (particles[i].sourceTerm - particles[i].Ap) / particles[i].diagElem);
          } else {
            particles[i].pressure = 0;
          }
        }
      }
      // Calculate density error considering only fluid particles
      avgDensError = 0.0;
      for (size_t i = 0; i < numParticles; i++) {
        if (particles[i].fluid) {
          avgDensError += std::abs(particles[i].Ap - particles[i].sourceTerm);
        }
      }
      avgDensError /= params.fluidHeight * params.fluidWidth;
      errorInPercent = avgDensError / params.restDensity;

      itJacobi++;
    }
    avgJacobiIterations += itJacobi;
    advectParticles(particles, numParticles, params.deltaT);

    auto endSolver = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> solverDuration = endSolver - startSolver;
    std::chrono::duration<double> total = neighborSearchDuration + initDuration + solverDuration;

    // Output outside time measurement
    std::cout << itJacobi << "  " << errorInPercent * 100 << "  " << params.toleratedError << std::endl;

    if (it % 100 == 0) {
      std::cout << "Avg number of Jacobi iterations: " << avgJacobiIterations / static_cast<double>(it) << std::endl;
    }

    totalSolverTime += solverDuration;

    if (it % params.snapshotFrequency == 0) {
      iisph::visualization::renderSnapshot(window, particles, numParticles, params.h, simHeight, it);
    }

    if (it == params.saveState) {
      iisph::io::saveParticlesToCSV(particles, numParticles, "particle_state");
    }
  }

  std::cout << "Gesamtlaufzeit Jacobi init: " << totalInitTime.count() << std::endl;
  std::cout << "Gesamtlaufzeit Jacobi solver: " << totalSolverTime.count() << std::endl;
  std::cout << "Gesamtlaufzeit Nachbarsuche: " << totalNeighborTime.count() << std::endl;
  std::cout << "Insgesamt: " << totalNeighborTime.count() + totalInitTime.count() + totalSolverTime.count()
            << std::endl;

  delete[] particles;
  delete[] neighbors;
  return 0;
}
