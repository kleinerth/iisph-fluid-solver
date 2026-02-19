#include "fluid/neighborSearch.h"
#include "fluid/particleBasics.h"
#include <algorithm> // find
#include <cmath>     // ceil, floor
#include <vector>    // vector, size_t

namespace iisph::fluid {

void countingSort(iisph::fluid::Particle* particles, const size_t numParticles, int digitPlace) {
  iisph::fluid::Particle* output = new iisph::fluid::Particle[numParticles];
  int count[10] = {0};

  for (size_t i = 0; i < numParticles; ++i) {
    int digit = (particles[i].cell_ID / digitPlace) % 10;
    count[digit]++;
  }

  for (int i = 1; i < 10; ++i) {
    count[i] += count[i - 1];
  }

  for (int i = static_cast<int>(numParticles) - 1; i >= 0; --i) {
    int digit = (particles[i].cell_ID / digitPlace) % 10;
    output[count[digit] - 1] = particles[i];
    count[digit]--;
  }

  for (size_t i = 0; i < numParticles; ++i) {
    particles[i] = output[i];
  }

  delete[] output;
}

void radixSort(iisph::fluid::Particle* particles, const size_t numParticles) {
  int max = particles[0].cell_ID;
  for (size_t i = 1; i < numParticles; ++i) {
    if (particles[i].cell_ID > max) {
      max = particles[i].cell_ID;
    }
  }
  for (int digitPlace = 1; max / digitPlace > 0; digitPlace *= 10) {
    countingSort(particles, numParticles, digitPlace);
  }
}

void insertionSort(iisph::fluid::Particle* particles, const size_t numParticles) {
  for (size_t i = 1; i < numParticles; ++i) {
    iisph::fluid::Particle key = particles[i];
    size_t j = i;
    while (j > 0 && particles[j - 1].cell_ID > key.cell_ID) {
      particles[j] = particles[j - 1];
      --j;
    }
    particles[j] = key;
  }
}

int binarySearch(const int* filledCells, const int arrSize, const int target) {
  int low = 0, high = arrSize - 1;
  while (low <= high) {
    int mid = (high - low) / 2 + low;
    if (filledCells[mid] == target) {
      return mid + 1;
    }
    if (filledCells[mid] > target) {
      high = mid - 1;
    } else {
      low = mid + 1;
    }
  }
  return 0;
}

void zIndexSort(iisph::fluid::Particle* particles, std::vector<int>* neighbors, size_t numParticles, double h,
                double xmin, double xmax, double ymin, double ymax, int sortingMethod) {

  int K = std::ceil((xmax - xmin) / (2 * h));

  // Calculate xyz-index
  for (size_t i = 0; i < numParticles; i++) {
    int k = static_cast<int>(floor((particles[i].posX) / (2 * h)));
    int l = static_cast<int>(floor((particles[i].posY) / (2 * h)));
    particles[i].cell_ID = k + l * K;
  }

  // Sort particles in-place with respect to their cell ID
  if (sortingMethod == 1) {
    insertionSort(particles, numParticles);
  }
  if (sortingMethod == 2) {
    radixSort(particles, numParticles);
  }

  // Determine size of compact cell array by counting the number of different
  // cell indices
  int numFilledCells = 1;
  for (size_t i = 0; i < numParticles - 1; i++) {
    if (particles[i].cell_ID != particles[i + 1].cell_ID) {
      numFilledCells++;
    }
  }

  // Create two info arrays. In order to calculate the number of particles in
  // the last cell, the second array needs one more entry than the paper
  // mentions
  int* firstParticleInCell = new int[numFilledCells + 1]();
  int* filledCells = new int[numFilledCells]();

  int marker = 1;
  int scan = 1;
  firstParticleInCell[0] = 0;
  filledCells[0] = particles[0].cell_ID;
  for (size_t i = 1; i < numParticles; i++) {
    if (particles[i].cell_ID != particles[i - 1].cell_ID) {
      marker = 1;
      firstParticleInCell[scan] = i;
      filledCells[scan] = particles[i].cell_ID;
    } else {
      marker = 0;
    }
    scan += marker;
  }
  firstParticleInCell[numFilledCells] = numParticles;

  // NEW QUERY (schwer lesbar - es gibt halt f�r eine range 3 F�lle:
  // 6 nicht gefunden, 7 nicht gefunden, 8 nicht gefunden.

  // Query
  // Is difficult to read, but minimizes the number of necessary searches.
  // Example: [6,7,8] cell range. Cases: 6 not found, but 7. 6,7 not found,
  // but 8.

  for (int c = 0; c < numFilledCells; c++) {
    int numParticlesInCell = firstParticleInCell[c + 1] - firstParticleInCell[c];

    if (numParticlesInCell > 0) {
      int first = firstParticleInCell[c];
      int k = static_cast<int>(floor((particles[first].posX) / (2 * h)));
      int l = static_cast<int>(floor((particles[first].posY) / (2 * h)));

      // Find neighbors of each particle in a certain cell
      for (int p = 0; p < numParticlesInCell; p++) {
        if (particles[firstParticleInCell[c] + p].fluid) {
          // Currently active particle: particles[firstParticleInCell[c] + p]
          // Clear the old neighbors of this particle
          neighbors[firstParticleInCell[c] + p].clear();
          // Iterate through all cells that neighbor the cell of the currently
          // active particle
          for (int y : {-1, 0, 1}) {
            int x = -1;
            int neighborCell = k + x + (l + y) * K;
            // Check if the currently active neighboring cell contains particles
            int found = binarySearch(filledCells, numFilledCells, neighborCell);
            if (found) {
              // 6 found
              int indexForPIC = found - 1; // where the current neighbor cell was found in filledCells
              int numParticlesNeighCell = firstParticleInCell[indexForPIC + 1] - firstParticleInCell[indexForPIC];
              for (int j = 0; j < numParticlesNeighCell; j++) {
                // Current possible neighbor candidate:
                // particles[firstParticleInCell[indexForPIC] + j]
                double dist = iisph::fluid::euclideanDistance(particles[firstParticleInCell[c] + p].posX,
                                                              particles[firstParticleInCell[c] + p].posY,
                                                              particles[firstParticleInCell[indexForPIC] + j].posX,
                                                              particles[firstParticleInCell[indexForPIC] + j].posY);
                if (dist < 2 * h) {
                  neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC] + j);
                }
              }
              // Current case: cell 6 was filled. Check if cell 7 also contains
              // particles by checking the next entry of filledCells
              if (filledCells[indexForPIC + 1] == neighborCell + 1) {
                // 7 found
                int numParticlesNeighCell = firstParticleInCell[indexForPIC + 2] - firstParticleInCell[indexForPIC + 1];
                for (int j = 0; j < numParticlesNeighCell; j++) {
                  // Current possible neighbor candidate:
                  // particles[firstParticleInCell[indexForPIC + 1] + j]
                  double dist = iisph::fluid::euclideanDistance(
                      particles[firstParticleInCell[c] + p].posX, particles[firstParticleInCell[c] + p].posY,
                      particles[firstParticleInCell[indexForPIC + 1] + j].posX,
                      particles[firstParticleInCell[indexForPIC + 1] + j].posY);
                  if (dist < 2 * h) {
                    neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC + 1] + j);
                  }
                }

                // Current case: cell 6 was filled, cell 7 was treated. Check if
                // cell 8 also contains particles.
                if (filledCells[indexForPIC + 2] == neighborCell + 2) {
                  int numParticlesNeighCell =
                      firstParticleInCell[indexForPIC + 3] - firstParticleInCell[indexForPIC + 2];
                  for (int j = 0; j < numParticlesNeighCell; j++) {
                    // Current possible neighbor candidate:
                    // particles[firstParticleInCell[indexForPIC + 2] + j]
                    double dist = iisph::fluid::euclideanDistance(
                        particles[firstParticleInCell[c] + p].posX, particles[firstParticleInCell[c] + p].posY,
                        particles[firstParticleInCell[indexForPIC + 2] + j].posX,
                        particles[firstParticleInCell[indexForPIC + 2] + j].posY);
                    if (dist < 2 * h) {
                      neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC + 2] + j);
                    }
                  }
                } else {
                  // 8 is empty
                }
              } else {
                // Current case: cell 6 was filled, 7 was empty. Check if cell 8
                // is next to cell 6, indicating that cell 8 is filled
                if (filledCells[indexForPIC + 1] == neighborCell + 2) {
                  int numParticlesNeighCell =
                      firstParticleInCell[indexForPIC + 2] - firstParticleInCell[indexForPIC + 1];
                  for (int j = 0; j < numParticlesNeighCell; j++) {
                    // Current possible neighbor candidate:
                    // particles[firstParticleInCell[indexForPIC + 1] + j]
                    double dist = iisph::fluid::euclideanDistance(
                        particles[firstParticleInCell[c] + p].posX, particles[firstParticleInCell[c] + p].posY,
                        particles[firstParticleInCell[indexForPIC + 1] + j].posX,
                        particles[firstParticleInCell[indexForPIC + 1] + j].posY);
                    if (dist < 2 * h) {
                      neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC + 1] + j);
                    }
                  }
                } else {
                  // 8 is empty
                }
              }

            } else {
              // 6 not found, ie cell 6 is empty. Search for 7
              found = binarySearch(filledCells, numFilledCells, neighborCell + 1);
              if (found) {
                int indexForPIC = found - 1;
                int numParticlesNeighCell = firstParticleInCell[indexForPIC + 1] - firstParticleInCell[indexForPIC];
                for (int j = 0; j < numParticlesNeighCell; j++) {
                  // Current possible neighbor candidate:
                  // particles[firstParticleInCell[indexForPIC] + j]
                  double dist = iisph::fluid::euclideanDistance(particles[firstParticleInCell[c] + p].posX,
                                                                particles[firstParticleInCell[c] + p].posY,
                                                                particles[firstParticleInCell[indexForPIC] + j].posX,
                                                                particles[firstParticleInCell[indexForPIC] + j].posY);
                  if (dist < 2 * h) {
                    neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC] + j);
                  }
                }
                // Check if cell 8 is next entry of filledCells
                if (filledCells[indexForPIC + 1] == neighborCell + 2) {
                  int numParticlesNeighCell =
                      firstParticleInCell[indexForPIC + 2] - firstParticleInCell[indexForPIC + 1];
                  for (int j = 0; j < numParticlesNeighCell; j++) {
                    // Current possible neighbor candidate:
                    // particles[firstParticleInCell[indexForPIC + 1] + j]
                    double dist = iisph::fluid::euclideanDistance(
                        particles[firstParticleInCell[c] + p].posX, particles[firstParticleInCell[c] + p].posY,
                        particles[firstParticleInCell[indexForPIC + 1] + j].posX,
                        particles[firstParticleInCell[indexForPIC + 1] + j].posY);
                    if (dist < 2 * h) {
                      neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC + 1] + j);
                    }
                  }
                } else {
                  // 8 empty
                }
              } else {
                // 6 and 7 empty. Search for 8
                found = binarySearch(filledCells, numFilledCells, neighborCell + 2);
                if (found) {
                  int indexForPIC = found - 1;
                  int numParticlesNeighCell = firstParticleInCell[indexForPIC + 1] - firstParticleInCell[indexForPIC];
                  for (int j = 0; j < numParticlesNeighCell; j++) {
                    double dist = iisph::fluid::euclideanDistance(particles[firstParticleInCell[c] + p].posX,
                                                                  particles[firstParticleInCell[c] + p].posY,
                                                                  particles[firstParticleInCell[indexForPIC] + j].posX,
                                                                  particles[firstParticleInCell[indexForPIC] + j].posY);
                    if (dist < 2 * h) {
                      neighbors[firstParticleInCell[c] + p].push_back(firstParticleInCell[indexForPIC] + j);
                    }
                  }
                } else {
                  // 8 empty
                }
              }
            }
          }
        }
      }
    }
  }

  delete[] firstParticleInCell;
  delete[] filledCells;
}

} // namespace iisph::fluid
