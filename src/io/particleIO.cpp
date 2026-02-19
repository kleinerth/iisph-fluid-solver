#include "io/particleIO.h"
#include "fluid/particleBasics.h"
#include <iostream>
#include <fstream> // std::ifstream, std::ofstream
#include <sstream> // std::stringstream
#include <string>

namespace iisph::io {

bool saveParticlesToCSV(iisph::fluid::Particle* particles, size_t numParticles, const std::string& filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: file did not open. The particle state was not saved.\n";
    return false;
  }

  file << "posX,posY,density,pressure,mass,vx,vy,vxPredicted,vyPredicted,ax,ay,"
          "sourceterm,diagelem,fluid,cell_ID\n";

  for (size_t i = 0; i < numParticles; ++i) {
    const iisph::fluid::Particle& p = particles[i];
    file << p.posX << "," << p.posY << "," << p.density << "," << p.pressure << "," << p.mass << "," << p.vx << ","
         << p.vy << "," << p.vxPredicted << "," << p.vyPredicted << "," << p.ax << "," << p.ay << "," << p.sourceTerm
         << "," << p.diagElem << "," << p.fluid << "," << p.cell_ID << "\n";
  }

  std::cout << "Particle state successfully saved" << std::endl;

  file.close();

  return true;
}

bool loadParticlesFromCSV(iisph::fluid::Particle*& particles, size_t numParticles, const std::string& filename) {
  particles = new iisph::fluid::Particle[numParticles];

  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: particle state could not be opened.\n";
    return false;
  }

  std::string line;
  std::getline(file, line); // Skip header

  for (size_t i = 0; i < numParticles; ++i) {
    if (!std::getline(file, line)) {
      std::cerr << "Error: number of particles in simParameters and saved "
                   "state does not match \n";
      delete[] particles;
      return false;
    }

    std::stringstream lineStream(line);
    std::string value;
    iisph::fluid::Particle& p = particles[i];

    std::getline(lineStream, value, ',');
    p.posX = std::stod(value);
    std::getline(lineStream, value, ',');
    p.posY = std::stod(value);
    std::getline(lineStream, value, ',');
    p.density = std::stod(value);
    std::getline(lineStream, value, ',');
    p.pressure = std::stod(value);
    std::getline(lineStream, value, ',');
    p.mass = std::stod(value);
    std::getline(lineStream, value, ',');
    p.vx = std::stod(value);
    std::getline(lineStream, value, ',');
    p.vy = std::stod(value);
    std::getline(lineStream, value, ',');
    p.vxPredicted = std::stod(value);
    std::getline(lineStream, value, ',');
    p.vyPredicted = std::stod(value);
    std::getline(lineStream, value, ',');
    p.ax = std::stod(value);
    std::getline(lineStream, value, ',');
    p.ay = std::stod(value);
    std::getline(lineStream, value, ',');
    p.sourceTerm = std::stod(value);
    std::getline(lineStream, value, ',');
    p.diagElem = std::stod(value);
    std::getline(lineStream, value, ',');
    p.fluid = std::stoi(value);
    std::getline(lineStream, value, ',');
    p.cell_ID = std::stoi(value);
  }

  std::cout << "Successfully loaded particle state" << std::endl;
  return true;
}

} // namespace iisph::io
