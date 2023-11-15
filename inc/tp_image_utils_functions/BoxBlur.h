#pragma once

#include <vector>
#include <cmath>
#include <thread>
#include <atomic>

#include <glm/glm.hpp>

namespace tp_image_utils_functions{

//##################################################################################################
std::vector<size_t> boxesForGauss(float sigma, size_t n);

//##################################################################################################
void gaussBlur_4_cpu(float* scl, float* aux, size_t w, size_t h, size_t r);

}

