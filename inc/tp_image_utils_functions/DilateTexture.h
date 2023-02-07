#ifndef tp_image_utils_functions_DilateTexture_h
#define tp_image_utils_functions_DilateTexture_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/ByteMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
bool TP_IMAGE_UTILS_FUNCTIONS_EXPORT dilateTexture(const tp_image_utils::ByteMap& mask,
                                                   tp_image_utils::ColorMap& texture,
                                                   size_t radius);

}

#endif
