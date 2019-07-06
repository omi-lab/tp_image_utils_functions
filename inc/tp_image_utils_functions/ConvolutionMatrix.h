#ifndef tp_image_utils_functions_ConvolutionMatrix_h
#define tp_image_utils_functions_ConvolutionMatrix_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/ByteMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
class ConvolutionMatrix
{
public:

  //################################################################################################
  ConvolutionMatrix();

  //################################################################################################
  ConvolutionMatrix(const std::string& text);

  //################################################################################################
  ConvolutionMatrix(const std::vector<double>& matrixData, size_t width, size_t height);

  //################################################################################################
  size_t width()const;

  //################################################################################################
  size_t height()const;

  //################################################################################################
  const std::vector<double>& matrixData()const;

  //################################################################################################
  void setMatrixData(const std::vector<double>& matrixData, size_t width, size_t height);

  //################################################################################################
  std::string toString()const;

  //################################################################################################
  void loadString(const std::string& text);

  //################################################################################################
  tp_image_utils::ColorMap convolve(const tp_image_utils::ColorMap& src)const;

  //################################################################################################
  //! Create a 3x3 identity matrix
  void makeIdentity();

  //################################################################################################
  void makeBlur();

private:
  std::vector<double> m_matrixData;
  size_t m_width{0};
  size_t m_height{0};
};

//##################################################################################################
//! Extract the colors and remove shading
/*!
Apply a convolution matrix to the image and return the result.

\note the width and height should be odd numbers larger than one.

\param src - The source image.
\param matrixData - The matrix organised as rows.
\param width - The number of columns in the matrix.
\param height - The number of rows in the matrix.

\return The image with the convolution matrix applied.
*/
tp_image_utils::ColorMap convolutionMatrix(const tp_image_utils::ColorMap& src, const std::vector<double>& matrixData, size_t width, size_t height);
}

#endif
