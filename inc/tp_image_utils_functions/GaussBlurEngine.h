#pragma once

#include <string>
#include <memory>

//#define USE_OPENCL

namespace tp_image_utils_functions
{

//##################################################################################################
class GaussBlurEngine{
public:
  GaussBlurEngine();
  ~GaussBlurEngine();

  void doBlur(float* scl, size_t w, size_t h, size_t r);

  std::string getErrorString();
  std::string getInfoString();

private:
  class GaussBlurAccelerator;
  std::unique_ptr<GaussBlurAccelerator> oclGaussBlur;
};

}
