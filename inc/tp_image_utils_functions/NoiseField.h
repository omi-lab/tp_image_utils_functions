#ifndef tp_image_utils_functions_NoiseField_h
#define tp_image_utils_functions_NoiseField_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/ByteMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ByteMap noiseField(const tp_image_utils::ByteMap& src, int radius);

//##################################################################################################
tp_image_utils::ByteMap noiseFieldGrid(const tp_image_utils::ByteMap& src, int cellSize);

}

#endif
