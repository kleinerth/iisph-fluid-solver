#include "visualization/graphics.h"
#include "fluid/particleBasics.h"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>    // sqrt
#include <sstream>
#include <string>

namespace iisph::visualization {

float setupWindow(sf::RenderWindow& window, const size_t baseLength, const double h, const size_t wallHeight,
                  const size_t wallThickness) {
  float simWidth = baseLength * h;
  float simHeight = (wallHeight + wallThickness + 2) * h + h;
  sf::FloatRect simBounds(0, 0, simWidth, simHeight);

  sf::View view;
  float windowRatio = static_cast<float>(window.getSize().x) / window.getSize().y;
  float simRatio = simBounds.width / simBounds.height;

  if (windowRatio > simRatio) {
    float newWidth = simBounds.height * windowRatio;
    view.setSize(newWidth, simBounds.height);
  } else {
    float newHeight = simBounds.width / windowRatio;
    view.setSize(simBounds.width, newHeight);
  }

  view.setCenter(simBounds.width / 2.f, simBounds.height / 2.f);
  window.setView(view);

  return simHeight;
}

void saveWindowAsImage(sf::RenderWindow& window, const std::string& filename) {
  sf::Image screenshot = window.capture();
  if (!screenshot.saveToFile(filename)) {
    std::cerr << "Error: Unable to save screenshot to " << filename << std::endl;
  }
}

void renderSnapshot(sf::RenderWindow& window, const iisph::fluid::Particle* particles, const size_t numParticles,
                    const double h, const float simHeight, const int it) {
  window.clear(sf::Color::White);

  float particleRadius = h / 2.0f;
  float maxVelocityForPlot = 8.0f;

  for (size_t i = 0; i < numParticles; i++) {
    sf::CircleShape circle(particleRadius);
    circle.setOrigin(particleRadius, particleRadius);
    circle.setPosition(particles[i].posX,
                       simHeight - particles[i].posY); // Flip y-axis

    // Color = velocity
    if (particles[i].fluid) {
      float velocity = std::sqrt(particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
      float normalizedVelocity = std::min(velocity / maxVelocityForPlot, 1.0f);
      int red = std::min(255, static_cast<int>(normalizedVelocity * 255));
      int blue = 255 - red;
      circle.setFillColor(sf::Color(red, 0, blue));
    } else {
      circle.setFillColor(sf::Color(127, 121, 121));
    }
    window.draw(circle);
  }
  window.display();

  std::stringstream filename;
  filename << "output_images/"
           << "img_" << it << ".jpg";
  saveWindowAsImage(window, filename.str());
}

} // namespace iisph::visualization
