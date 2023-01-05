#ifndef tp_image_utils_functions_ConvolutionMatrix_h
#define tp_image_utils_functions_ConvolutionMatrix_h

#include "tp_image_utils_functions/Globals.h" // IWYU pragma: keep

#include "glm/glm.hpp" // IWYU pragma: keep

namespace tp_image_utils
{
class ColorMap;
class ColorMapF;
}

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
  [[nodiscard]]size_t width()const;

  //################################################################################################
  [[nodiscard]]size_t height()const;

  //################################################################################################
  const std::vector<double>& matrixData()const;

  //################################################################################################
  std::vector<float> matrixDataF()const;

  //################################################################################################
  void setMatrixData(const std::vector<double>& matrixData, size_t width, size_t height);

  //################################################################################################
  [[nodiscard]]std::string toString()const;

  //################################################################################################
  void loadString(const std::string& text);

  //################################################################################################
  [[nodiscard]]tp_image_utils::ColorMap convolve(const tp_image_utils::ColorMap& src)const;

  //################################################################################################
  //! Create a 3x3 identity matrix
  void makeIdentity();

  //################################################################################################
  void makeBlur3();

  //################################################################################################
  void makeBlur5();

  //################################################################################################
  void divideBySize();

private:
  std::vector<double> m_matrixData;
  size_t m_width{0};
  size_t m_height{0};
};

//##################################################################################################
//! Apply a convolution matrix to the image
/*!
Apply a convolution matrix to the image and return the result.

\note the width and height should be odd numbers larger than one.

\param src - The source image.
\param matrixData - The matrix organised as rows.
\param width - The number of columns in the matrix.
\param height - The number of rows in the matrix.

\return The image with the convolution matrix applied.
*/
tp_image_utils::ColorMap convolutionMatrix(const tp_image_utils::ColorMap& src,
                                           const std::vector<double>& matrixData,
                                           size_t width,
                                           size_t height);

//##################################################################################################
tp_image_utils::ColorMapF convolvePadded(const tp_image_utils::ColorMapF& src,
                                         const std::vector<float>& matrixData,
                                         size_t width,
                                         size_t height);

//##################################################################################################
std::vector<glm::vec3> gaussBlur(std::vector<glm::vec3>& source,
                                 size_t w,
                                 size_t h,
                                 size_t radius);

//##################################################################################################
tp_image_utils::ColorMap gaussBlur(const tp_image_utils::ColorMap& source, size_t radius);

//##################################################################################################
tp_image_utils::ColorMapF gaussBlur(const tp_image_utils::ColorMapF& source, size_t radius);

//##################################################################################################
tp_image_utils::ColorMapF blur3(const tp_image_utils::ColorMapF& src);

//##################################################################################################
tp_image_utils::ColorMapF blur5(const tp_image_utils::ColorMapF& src);
}

#endif
