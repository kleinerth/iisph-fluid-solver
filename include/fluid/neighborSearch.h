#ifndef NEIGHBORSEARCH_H
#define NEIGHBORSEARCH_H

#include "fluid/particleBasics.h"
#include <vector>    // vector, size_t

namespace iisph::fluid {

void countingSort(iisph::fluid::Particle* particles, const size_t numParticles, int digitPlace);
void radixSort(iisph::fluid::Particle* particles, const size_t numParticles);
void insertionSort(iisph::fluid::Particle* particles, const size_t numParticles);
int binarySearch(const int* filledCells, const int arrSize, const int target);
void zIndexSort(iisph::fluid::Particle* particles, std::vector<int>* neighbors, size_t numParticles, double h,
                double xmin, double xmax, double ymin, double ymax, int sortingMethod);

} // namespace iisph::fluid

#endif
