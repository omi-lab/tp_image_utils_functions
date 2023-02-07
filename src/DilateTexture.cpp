#include "tp_image_utils_functions/DilateTexture.h"

#include "tp_image_utils/ColorMap.h"

#include "glm/glm.hpp" // IWYU pragma: keep

namespace tp_image_utils_functions
{

//##################################################################################################
bool dilateTexture(const tp_image_utils::ByteMap& mask,
                   tp_image_utils::ColorMap& texture,
                   size_t radius)
{
  size_t w = texture.width ();
  size_t h = texture.height();

  if(mask.width() != w || mask.height() != h)
    return false;

  tp_image_utils::ByteMap mask1 = mask;
  tp_image_utils::ByteMap mask2 = mask;
  for(size_t pass=0; pass<radius; pass++)
  {
    auto m = mask.constData();
    auto n = mask2.data();
    auto p = texture.data();

    for(size_t y=0; y<h; y++)
    {
      for(size_t x=0; x<w; x++, m++, n++, p++)
      {
        if((*m) == 0)
          continue;

        auto tryPixel = [&](ptrdiff_t xo, ptrdiff_t yo)
        {
          size_t xoo = size_t(ptrdiff_t(x)-xo);
          size_t yoo = size_t(ptrdiff_t(y)-yo);

          if(mask1.pixel(xoo, yoo, 255) == 0)
          {
            (*p) = texture.pixel(xoo, yoo);
            (*n) = 0;
            return true;
          }
          return false;
        };

        using V = glm::vec<2, ptrdiff_t, glm::defaultp>;
        for(auto i : {V(-1,0), V(0,-1), V(1,0), V(0,1), V(-1,-1), V(1,-1), V(1,1), V(-1,1)})
          if(tryPixel(i.x, i.y))
            break;
      }
    }

    mask1 = mask2;
  }

  return true;
}

}
