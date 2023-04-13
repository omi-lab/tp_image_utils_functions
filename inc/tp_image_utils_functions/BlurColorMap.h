#ifndef tp_image_utils_functions_BlurColorMap_h
#define tp_image_utils_functions_BlurColorMap_h

#include "tp_image_utils/ColorMapF.h"


namespace tp_image_utils_functions
{
void blurColorMap(tp_image_utils::ColorMapF& colorMap, size_t radius);
}

#endif
