#include "tp_image_utils_functions/BlurColorMap.h"

#include "tp_image_utils_functions/ConvolutionMatrix.h"

#include "tp_utils/Parallel.h"

#include <cstring>

namespace tp_image_utils_functions
{

//##################################################################################################
void blurColorMap(tp_image_utils::ColorMapF& colorMap, size_t radius)
{
  if(radius<1)
    return;

  size_t width  = colorMap.width ();
  size_t height = colorMap.height();

  std::vector<glm::vec3> input((width*height));

  size_t xOffset = colorMap.width ();
  size_t yOffset = colorMap.height();

  {
    size_t c=0;
    tp_utils::parallel([&](auto locker)
    {
      for(;;)
      {
        size_t y;
        locker([&]{y=c; c++;});

        if(y>=height)
          return;

        size_t ySource = (y+yOffset) % colorMap.height();
        for(size_t x=0; x<width; x++)
        {
          size_t xSource = (x+xOffset) % colorMap.width();
          input[x+(y*width)] = glm::vec4(colorMap.pixel(xSource, ySource));
        }
      }
    });
  }

  std::vector<glm::vec3> output = tp_image_utils_functions::gaussBlur(input, colorMap.width(), colorMap.height(), radius);

  {
    glm::vec4* dstData = colorMap.data();

    size_t c=0;
    tp_utils::parallel([&](auto locker)
    {
      for(;;)
      {
        size_t y;
        locker([&]{y=c; c++;});

        if(y>=colorMap.height())
          return;

        glm::vec4* dst = dstData + (y*colorMap.width());
        glm::vec4* dstMax = dst+colorMap.width();
        size_t offset = ((y) * colorMap.width() );
        for(; dst<dstMax; dst++, offset++)
          *dst = glm::vec4(output[offset], 1.0f);
      }
    });
  }
}

}
