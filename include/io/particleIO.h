#ifndef PARTICLEIO_H
#define PARTICLEIO_H

#include "fluid/particleBasics.h"
#include <string>

namespace iisph::io {

bool saveParticlesToCSV(iisph::fluid::Particle* particles, size_t numParticles, const std::string& filename);
bool loadParticlesFromCSV(iisph::fluid::Particle*& particles, size_t numParticles, const std::string& filename);

} // namespace iisph::io

#endif
