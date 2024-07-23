#ifndef tp_image_utils_functions_DiffImages_h
#define tp_image_utils_functions_DiffImages_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/ColorMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
bool TP_IMAGE_UTILS_FUNCTIONS_EXPORT diffImages(const tp_image_utils::ColorMap& a,
                                                const tp_image_utils::ColorMap& b,
                                                double amplification,
                                                tp_image_utils::ColorMap& result);
}

#endif
