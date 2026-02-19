#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "fluid/particleBasics.h"
#include <SFML/Graphics.hpp>
#include <string>

namespace iisph::visualization {

float setupWindow(sf::RenderWindow& window, const size_t baseLength, const double h, const size_t wallHeight,
                  const size_t wallThickness);
void saveWindowAsImage(sf::RenderWindow& window, const std::string& filename);
void renderSnapshot(sf::RenderWindow& window, const iisph::fluid::Particle* particles, const size_t numParticles,
                    const double h, const float simHeight, const int it);

} // namespace iisph::visualization

#endif
