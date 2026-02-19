#include "fluid/simMethods.h"
#include "fluid/kernel.h"
#include "fluid/particleBasics.h"
#include <vector>

namespace iisph::fluid {

void computeDensity(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
                    const double h, const double alpha, const double gammaDens) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      double dens = 0;
      for (int neighbor : neighbors[i]) {
        if (particles[neighbor].fluid) {
          dens += particles[neighbor].mass *
                  W(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        } else {
          dens += gammaDens * particles[neighbor].mass *
                  W(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        }
      }
      particles[i].density = dens;
    }
  }
}

void computeNonPressureAcc(iisph::fluid::Particle* particles, const size_t numParticles,
                           const std::vector<int>* neighbors, const double h, const double alpha, const double ny,
                           const double gravity) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      double viscX = 0;
      double viscY = 0;
      for (int neighbor : neighbors[i]) {
        std::pair<double, double> gradient =
            gradW(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        double dx = particles[i].posX - particles[neighbor].posX;
        double dy = particles[i].posY - particles[neighbor].posY;
        double dvx = particles[i].vx - particles[neighbor].vx;
        double dvy = particles[i].vy - particles[neighbor].vy;
        double factor = 2 * ny * (particles[neighbor].mass / particles[neighbor].density) *
                        ((dvx * dx + dvy * dy) / (dx * dx + dy * dy + 0.01 * h * h));
        viscX += factor * gradient.first;
        viscY += factor * gradient.second;
      }
      viscY += gravity;
      particles[i].ax = viscX;
      particles[i].ay = viscY;
    }
  }
}

void predictVelocity(iisph::fluid::Particle* particles, const size_t numParticles, const double deltaT) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      particles[i].vxPredicted = particles[i].vx + deltaT * particles[i].ax;
      particles[i].vyPredicted = particles[i].vy + deltaT * particles[i].ay;
    }
  }
}

void computeSourceTerm(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
                       const double h, const double alpha, const double restDensity, const double deltaT) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      double sourcet = 0;
      for (int neighbor : neighbors[i]) {
        std::pair<double, double> gradient =
            gradW(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        if (particles[neighbor].fluid) {
          sourcet += deltaT * particles[neighbor].mass *
                     ((particles[i].vxPredicted - particles[neighbor].vxPredicted) * gradient.first +
                      (particles[i].vyPredicted - particles[neighbor].vyPredicted) * gradient.second);
        } else {
          sourcet += deltaT * particles[neighbor].mass *
                     (particles[i].vxPredicted * gradient.first + particles[i].vyPredicted * gradient.second);
        }
      }
      particles[i].sourceTerm = restDensity - particles[i].density - sourcet;
    }
  }
}

void computeDiagonalElement(iisph::fluid::Particle* particles, const size_t numParticles,
                            const std::vector<int>* neighbors, const double h, const double alpha, const double gamma,
                            const double deltaT) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      // c_i
      double cx = 0;
      double cy = 0;
      for (int neighbor : neighbors[i]) {
        std::pair<double, double> gradient =
            gradW(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        if (particles[neighbor].fluid) {
          cx -= particles[neighbor].mass / (particles[i].density * particles[i].density) * gradient.first;
          cy -= particles[neighbor].mass / (particles[i].density * particles[i].density) * gradient.second;
        } else {
          cx -= 2 * gamma * particles[neighbor].mass / (particles[i].density * particles[i].density) * gradient.first;
          cy -= 2 * gamma * particles[neighbor].mass / (particles[i].density * particles[i].density) * gradient.second;
        }
      }

      double aii = 0;
      for (int neighbor : neighbors[i]) {
        std::pair<double, double> gradient =
            gradW(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        std::pair<double, double> gradient2 =
            gradW(particles[neighbor].posX, particles[neighbor].posY, particles[i].posX, particles[i].posY, h, alpha);
        if (particles[neighbor].fluid) {
          aii += particles[neighbor].mass * (cx * gradient.first + cy * gradient.second) +
                 particles[neighbor].mass * particles[i].mass / (particles[i].density * particles[i].density) *
                     ((gradient2.first * gradient.first) + (gradient2.second * gradient.second));
        } else {
          aii += particles[neighbor].mass * (cx * gradient.first + cy * gradient.second);
        }
      }
      particles[i].diagElem = deltaT * deltaT * aii;
    }
  }
}

void computePressureAcc(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
                        const double h, const double alpha, const double gamma) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      double pressAccX = 0.0;
      double pressAccY = 0.0;
      for (int neighbor : neighbors[i]) {
        std::pair<double, double> gradient =
            gradW(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        if (particles[neighbor].fluid) {
          double factor = particles[neighbor].mass *
                          (particles[i].pressure / (particles[i].density * particles[i].density) +
                           particles[neighbor].pressure / (particles[neighbor].density * particles[neighbor].density));
          pressAccX -= factor * gradient.first;
          pressAccY -= factor * gradient.second;
        } else {
          double factor =
              particles[neighbor].mass * 2 * particles[i].pressure / (particles[i].density * particles[i].density);
          pressAccX -= gamma * factor * gradient.first;
          pressAccY -= gamma * factor * gradient.second;
        }
      }
      particles[i].ax = pressAccX;
      particles[i].ay = pressAccY;
    }
  }
}

void computeAp(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
               const double h, double const alpha, const double deltaT) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      double productAp = 0;
      for (int neighbor : neighbors[i]) {
        std::pair<double, double> gradient =
            gradW(particles[i].posX, particles[i].posY, particles[neighbor].posX, particles[neighbor].posY, h, alpha);
        if (particles[neighbor].fluid) {
          productAp += deltaT * deltaT * particles[neighbor].mass *
                       ((particles[i].ax - particles[neighbor].ax) * gradient.first +
                        (particles[i].ay - particles[neighbor].ay) * gradient.second);
        } else {
          productAp += deltaT * deltaT * particles[neighbor].mass *
                       (particles[i].ax * gradient.first + particles[i].ay * gradient.second);
        }
      }
      particles[i].Ap = productAp;
    }
  }
}

void advectParticles(iisph::fluid::Particle* particles, const size_t numParticles, const double deltaT) {
  for (size_t i = 0; i < numParticles; i++) {
    if (particles[i].fluid) {
      particles[i].vx = particles[i].vxPredicted + deltaT * particles[i].ax;
      particles[i].vy = particles[i].vyPredicted + deltaT * particles[i].ay;
      particles[i].posX += deltaT * particles[i].vx;
      particles[i].posY += deltaT * particles[i].vy;
    }
  }
}

} // namespace iisph::fluid
