#include "tp_image_utils_functions/RotateSkyBox.h"

#include "tp_utils/Globals.h"

#include "glm/gtc/constants.hpp" // for pi

#include "tinyexr.h"

namespace tp_image_utils_functions
{

//##################################################################################################
bool rotateSkyBox(const std::string& inputHDRIPath, const glm::mat3& R, std::string& outputHDRIPath, std::string& errorMessage)
{
  int width=0;
  int height=0;
  float* rgba=nullptr;
  char* err=nullptr;
  TP_CLEANUP([&]
  {
    free(rgba);
    free(err);
  });

  if(LoadEXR(&rgba, &width, &height, inputHDRIPath.c_str(), const_cast<const char **>(&err)) != 0)
  {
    errorMessage = "Failed to load HDRI file " + inputHDRIPath + ": " + err;
    return false;
  }

  if(width<0 || height<0 || width>32768 || height>32768)
  {
    errorMessage = "Failed to load HDRI file " + inputHDRIPath + ": invalid image size";
    return false;
  }

  size_t w = size_t(width);
  size_t h = size_t(height);
  std::vector<float> outputRGBA(4*w*h);

  float invW = 1.f/float(w);
  float invH = 1.f/float(h);

  auto dst = outputRGBA.data();
  for(size_t y=0; y<h; y++)
    for(size_t x=0; x<w; x++, dst += 4)
    {
      glm::vec2 dstCoord{(0.5f+float(x))*invW, (0.5f+float(y))*invH};
      dstCoord.y = 1.f-dstCoord.y; // so that image coordinate y points up
      dstCoord -= 0.5f;
      dstCoord *= 2.0f; // x and y are now in range [-1,1]
      //dstCoord.x += 1.0f; // x is now in range [0,2], y range is [-1,1]

      float lon = dstCoord.x * glm::pi<float>(); // range [-pi,pi]
      float lat = dstCoord.y * (0.5f*glm::pi<float>()); // range [-pi/2,pi/2]

      glm::vec3 vec;
      vec.x = std::cos(lat) * std::sin(lon);
      vec.y = std::cos(lat) * std::cos(lon);
      vec.z = std::sin(lat);

      vec = R*vec;

      lon = std::atan2(vec.x, vec.y);
      lat = std::asin(vec.z);

      glm::vec2 srcCoord{lon/glm::pi<float>(), lat/(0.5f*glm::pi<float>())};
      //srcCoord.x -= 1.f;
      srcCoord *= 0.5f;
      srcCoord += 0.5f;
      srcCoord.y = 1.f-srcCoord.y;

      size_t srcX = std::clamp(int((srcCoord.x*float(w))+0.5f), 0, int(w-1));
      size_t srcY = std::clamp(int((srcCoord.y*float(h))+0.5f), 0, int(h-1));

      const float* srcPix = &rgba[4*w*srcY+4*srcX];
      for(int i=0; i<4; ++i)
        dst[i] = srcPix[i];
    }

  if(SaveEXR(outputRGBA.data(), width, height, 4/*RGBA*/, 0/*save_as_fp16*/, outputHDRIPath.c_str(), const_cast<const char **>(&err)) != 0)
  {
    errorMessage = "Failed to save HDRI file " + outputHDRIPath + ": " + err;
    return false;
  }

  return true;
}

}
