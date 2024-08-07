#ifndef tp_image_utils_functions_RotateSkyBox_h
#define tp_image_utils_functions_RotateSkyBox_h

#include "tp_image_utils_functions/Globals.h"

#include "glm/mat3x3.hpp"

#include <string>

namespace tp_image_utils_functions
{

//##################################################################################################
bool rotateSkyBox(const std::string& inputHDRIPath, const glm::mat3& R, std::string& outputHDRIPath, std::string& errorMessage);
  
}

#endif
