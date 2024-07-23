#include "tp_image_utils_functions/DiffImages.h"

namespace tp_image_utils_functions
{

//##################################################################################################
bool  diffImages(const tp_image_utils::ColorMap& a,
                 const tp_image_utils::ColorMap& b,
                 double amplification,
                 tp_image_utils::ColorMap& result)
{
  if(a.width() != b.width() || a.height() != b.height())
    return false;

  result.setSize(a.width(), a.height());

  auto srcA = a.constData();
  auto srcB = b.constData();
  auto dst = result.data();
  auto dstMax = dst + result.size();

  for(; dst<dstMax; srcA++, srcB++, dst++)
  {
    dst->r = uint8_t(std::clamp(amplification*std::fabs(int(srcA->r)-int(srcB->r)), 0., 255.));
    dst->g = uint8_t(std::clamp(amplification*std::fabs(int(srcA->g)-int(srcB->g)), 0., 255.));
    dst->b = uint8_t(std::clamp(amplification*std::fabs(int(srcA->b)-int(srcB->b)), 0., 255.));
    dst->a = 255;
  }

  return true;
}

}
