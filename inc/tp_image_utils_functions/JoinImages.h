#ifndef tp_image_utils_functions_JoinImages_h
#define tp_image_utils_functions_JoinImages_h

#include "tp_image_utils_functions/Globals.h" // IWYU pragma: keep

#include "tp_image_utils/ColorMap.h"

namespace tp_image_utils_functions
{
//##################################################################################################
tp_image_utils::ColorMap join2Images(const tp_image_utils::ColorMap& image1,
                                     const tp_image_utils::ColorMap& image2,
                                     size_t width,
                                     size_t height);

}

#endif
