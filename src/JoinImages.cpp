#include "tp_image_utils_functions/JoinImages.h"

#include <string.h>

namespace tp_image_utils_functions
{
//##################################################################################################
tp_image_utils::ColorMap join2Images(const tp_image_utils::ColorMap& image1,
                                     const tp_image_utils::ColorMap& image2,
                                     size_t width,
                                     size_t height)
{
  if(image1.size() != (width*height) ||
     image2.size() != (width*height))
    return tp_image_utils::ColorMap();

  tp_image_utils::ColorMap colorMap(width*2, height);

  auto s1 = image1.constData();
  auto s2 = image2.constData();
  auto dst = colorMap.data();
  for(size_t y=0; y<height; y++)
  {
    memcpy(dst, s1, width*4);
    dst+=width;
    s1+=width;
    memcpy(dst, s2, width*4);
    dst+=width;
    s2+=width;
  }

  return colorMap;
}

}
