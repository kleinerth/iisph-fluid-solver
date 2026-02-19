#ifndef SIMMETHODS_H
#define SIMMETHODS_H

#include "fluid/kernel.h"
#include "fluid/particleBasics.h"
#include <vector>

namespace iisph::fluid {

void computeDensity(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
                    const double h, const double alpha, const double gammaDens);
void computeNonPressureAcc(iisph::fluid::Particle* particles, const size_t numParticles,
                           const std::vector<int>* neighbors, const double h, const double alpha, const double ny,
                           const double gravity);
void predictVelocity(iisph::fluid::Particle* particles, const size_t numParticles, const double deltaT);
void computeSourceTerm(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
                       const double h, const double alpha, const double restDensity, const double deltaT);
void computeDiagonalElement(iisph::fluid::Particle* particles, const size_t numParticles,
                            const std::vector<int>* neighbors, const double h, const double alpha, const double gamma,
                            const double deltaT);
void computePressureAcc(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
                        const double h, const double alpha, const double gamma);
void computeAp(iisph::fluid::Particle* particles, const size_t numParticles, const std::vector<int>* neighbors,
               const double h, const double alpha, const double deltaT);
void advectParticles(iisph::fluid::Particle* particles, const size_t numParticles, const double deltaT);

} // namespace iisph::fluid

#endif
