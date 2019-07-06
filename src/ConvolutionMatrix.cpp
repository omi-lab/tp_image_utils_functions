#include "tp_image_utils_functions/ConvolutionMatrix.h"

#include "tp_image_utils/ColorMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
ConvolutionMatrix::ConvolutionMatrix()
{
  makeIdentity();
}

//##################################################################################################
ConvolutionMatrix::ConvolutionMatrix(const std::string& text)
{
  loadString(text);
}

//##################################################################################################
ConvolutionMatrix::ConvolutionMatrix(const std::vector<double>& matrixData, size_t width, size_t height)
{
  setMatrixData(matrixData, width, height);
}

//##################################################################################################
size_t ConvolutionMatrix::width()const
{
  return m_width;
}

//##################################################################################################
size_t ConvolutionMatrix::height()const
{
  return m_height;
}

//##################################################################################################
const std::vector<double>& ConvolutionMatrix::matrixData()const
{
  return m_matrixData;
}

//##################################################################################################
void ConvolutionMatrix::setMatrixData(const std::vector<double>& matrixData, size_t width, size_t height)
{
  if(width<3 || (width&1)!=1   ||
     height<3 || (height&1)!=1 ||
     (width*height) != matrixData.size())
    makeIdentity();
  else
  {
    m_width = width;
    m_height = height;
    m_matrixData = matrixData;
  }
}

//##################################################################################################
std::string ConvolutionMatrix::toString()const
{
  std::string output;

  const double* d = m_matrixData.data();

  for(size_t y=0; y<m_height; y++)
  {
    if(y>0)
      output.push_back('|');

    for(size_t x=0; x<m_width; x++)
    {
      if(x>0)
        output.push_back(',');

      output.append(std::to_string(*d));
      d++;
    }
  }

  return output;
}

//##################################################################################################
void ConvolutionMatrix::loadString(const std::string& text)
{
  m_matrixData.clear();
  m_height = 0;

  std::vector<std::string> rows;
  tpSplit(rows, text, '|', tp_utils::SplitBehavior::SkipEmptyParts);

  m_width = rows.size();

  if(m_width<3 || (m_width&1)!=1)
  {
    makeIdentity();
    return;
  }

  for(const std::string& row : tpConst(rows))
  {
    std::vector<std::string> parts;
    tpSplit(parts, row, ',', tp_utils::SplitBehavior::SkipEmptyParts);

    if(m_height<1)
    {
      m_height = parts.size();

      if(m_height<3 || (m_height&1)!=1)
      {
        makeIdentity();
        return;
      }
    }
    else if(m_height != parts.size())
    {
      makeIdentity();
      return;
    }

    try
    {
      for(const std::string& part : tpConst(parts))
        m_matrixData.push_back(std::stoi(part));
    }
    catch(...)
    {

    }
  }
}

//##################################################################################################
tp_image_utils::ColorMap ConvolutionMatrix::convolve(const tp_image_utils::ColorMap& src)const
{
  return convolutionMatrix(src, m_matrixData, m_width, m_height);
}

//##################################################################################################
void ConvolutionMatrix::makeIdentity()
{
  m_width =  3;
  m_height = 3;

  m_matrixData =
  {
    0, 0, 0,
    0, 9, 0,
    0, 0, 0
  };
}

//##################################################################################################
void ConvolutionMatrix::makeBlur()
{
  m_width =  5;
  m_height = 5;

  m_matrixData =
  {
    0, 0, 1, 0, 0,
    0, 1, 3, 1, 0,
    1, 3, 5, 3, 1,
    0, 1, 3, 1, 0,
    0, 0, 1, 0, 0
  };
}

namespace
{
//##################################################################################################
struct PixelDetails_lt
{
  double red  {0.0};
  double green{0.0};
  double blue {0.0};
};
}

//##################################################################################################
tp_image_utils::ColorMap convolutionMatrix(const tp_image_utils::ColorMap& src, const std::vector<double>& matrixData, size_t width, size_t height)
{
  if(width<=1 || height<=1)
    return tp_image_utils::ColorMap();

  size_t size = width*height;

  if(size>matrixData.size())
    return tp_image_utils::ColorMap();

  size_t dw = src.width()  - (width-1);//((width+1)/2);
  size_t dh = src.height() - (height-1);//((height+1)/2);

  if(dw<1 || dh<1 || dw>src.width() || dh>src.height())
    return tp_image_utils::ColorMap();

  std::vector<PixelDetails_lt> buffer;
  buffer.resize(dw*dh);

  PixelDetails_lt* bufferData = buffer.data();

  for(size_t my=0; my<height; my++)
  {
    for(size_t mx=0; mx<width;  mx++)
    {
      double weight = matrixData.at((my*width)+mx);
      for(size_t dy=0; dy<dh; dy++)
      {
        auto* s = src.constData() + ((dy+my)*src.width()) + mx;
        PixelDetails_lt* d = bufferData + (dy*dw);
        PixelDetails_lt* dMax = d + dw;
        while(d<dMax)
        {
          d->red   += double(s->r) * weight;
          d->green += double(s->g) * weight;
          d->blue  += double(s->b) * weight;

          d++;
          s++;
        }
      }
    }
  }

  tp_image_utils::ColorMap dst(dw, dh);
  {
    PixelDetails_lt* s = bufferData;

    for(size_t y=0; y<dh; y++)
    {
      TPPixel* d = dst.data() + (y*dw);
      TPPixel* dMax = d+dw;

      while(d<dMax)
      {
        d->r = uint8_t(tpBound(0, int(s->red  ) / int(size), 255));
        d->g = uint8_t(tpBound(0, int(s->green) / int(size), 255));
        d->b = uint8_t(tpBound(0, int(s->blue ) / int(size), 255));
        d->a = 255;
        d++;
        s++;
      }
    }
  }

  return dst;
}
}
