#include "tp_image_utils_functions/ToPolar.h"

#include <cmath>

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ByteMap toPolar(const tp_image_utils::ByteMap& src, size_t w, size_t h)
{
  if(w<2 || w>1024 || h<2 || h>1024)
    return tp_image_utils::ByteMap();

  if(src.width()<2 || src.height()<2)
    return tp_image_utils::ByteMap();

  tp_image_utils::ByteMap result(w, h);

  float fa = 6.28319f / float(h);
  float fx = 1.0f / float(w);
  float ox = float(src.width())  / 2.0f;
  float oy = float(src.height()) / 2.0f;

  for(size_t y=0; y<h; y++)
  {
    float a=float(y)*fa;
    float as = std::sin(a) * ox;
    float ac = std::cos(a) * oy;

    for(size_t x=0; x<w; x++)
    {
      float f = fx * float(x);
      auto sx = size_t(ox + (as*f));
      auto sy = size_t(oy + (ac*f));
      result.setPixel(x, y, src.pixel(sx, sy));
    }
  }

  return result;
}

}
