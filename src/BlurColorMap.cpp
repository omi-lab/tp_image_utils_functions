#include "tp_image_utils_functions/BlurColorMap.h"

#include "tp_image_utils_functions/ConvolutionMatrix.h"

#include "tp_utils/Parallel.h"

#include <cstring>
#include <atomic>

namespace tp_image_utils_functions
{

//##################################################################################################
void blurColorMap(tp_image_utils::ColorMapF& colorMap, size_t radius)
{
  if(radius<1)
    return;

  size_t const width  = colorMap.width ();
  size_t const height = colorMap.height();

  std::unique_ptr<glm::vec3[]> input( new glm::vec3[width*height]);
  std::unique_ptr<glm::vec3[]> output( new glm::vec3[width*height]);

  size_t xOffset = colorMap.width ();
  size_t yOffset = colorMap.height();


  {
    std::atomic<size_t> c{0};
    auto input_ptr = input.get();
    tp_utils::parallel([&](auto /*locker*/)
    {
      for(;;)
      {
        size_t y = c++;

        if(y>=height)
          return;

        size_t ySource = (y+yOffset) % colorMap.height();
        for(size_t x=0; x<width; x++)
        {
          size_t xSource = (x+xOffset) % colorMap.width();
          input_ptr[x+(y*width)] = colorMap.pixel(xSource, ySource);
        }
      }
    });
  }

  tp_image_utils_functions::gaussBlur(input.get(), output.get(), width, height, radius);

  {
    glm::vec4* dstData = colorMap.data();

    std::atomic<size_t> c{0};
    auto output_ptr = output.get();
    tp_utils::parallel([&](auto /*locker*/)
    {
      for(;;)
      {
        size_t y = c++;

        if(y>=colorMap.height())
          return;

        glm::vec4* dst = dstData + (y*colorMap.width());
        glm::vec4* dstMax = dst+colorMap.width();
        size_t offset = ((y) * colorMap.width() );
        for(; dst<dstMax; dst++, offset++)
          *dst = glm::vec4(output_ptr[offset], 1.0f);
      }
    });
  }
}

}
