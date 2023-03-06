#include "tp_image_utils_functions/ConvolutionMatrix.h"

#include "tp_image_utils/ColorMap.h"
#include "tp_image_utils/ColorMapF.h"

#include "tp_utils/Parallel.h"

#include <cstring>

namespace tp_image_utils_functions
{

namespace
{
//##################################################################################################
/*!
 * \param sigma standard deviation
 * \param n number of boxes
 * \return
 */
std::vector<size_t> boxesForGauss(float sigma, size_t n)
{
  // Ideal averaging filter width
  float wIdeal = std::sqrt((12.0f * sigma*sigma / float(n)) + 1.0f);

  size_t wl = size_t(std::floor(wIdeal));
  if(wl % 2 == 0)
    wl--;

  size_t wu = wl + 2;

  auto mIdeal = (12 * sigma*sigma - n * wl*wl - 4 * n*wl - 3 * n) / (-4 * wl - 4);
  size_t m = size_t(std::round(mIdeal));

  std::vector<size_t> sizes(n);
  for (size_t i=0; i<n; i++)
    sizes[i] = i < m ? wl : wu;

  return sizes;
}

//##################################################################################################
void boxBlurH_4(std::vector<glm::vec3>& scl, std::vector<glm::vec3>& tcl, size_t w, size_t h, size_t r)
{
  const size_t nHead = r+1;
  const size_t nTail = r;
  const size_t nBody = w-(nHead+nTail);

  const float iarr = 1.0f / float(r + r + 1);

  size_t c=0;
  tp_utils::parallel([&](auto locker)
  {
    for(;;)
    {
      size_t i;
      locker([&]{i=c; c++;});

      if(i>=h)
        return;

      {
        size_t ti = i * w;
        size_t li = ti;
        size_t ri = ti + r;

        glm::vec3* liSCL = &scl[li];
        glm::vec3* riSCL = &scl[ri];
        glm::vec3* tiTCL = &tcl[ti];

        glm::vec3 fv = scl[ti];
        glm::vec3 lv = scl[ti + w - 1];
        glm::vec3 val = float(r + 1)*fv;

        for(size_t j=0; j<r; j++)
          val += scl[ti+j];

        for(glm::vec3* tiTCLMax=tiTCL+nHead; tiTCL<tiTCLMax; riSCL++, tiTCL++)
        {
          val += *riSCL - fv;
          (*tiTCL) = val*iarr;
        }

        for(glm::vec3* tiTCLMax=tiTCL+nBody; tiTCL<tiTCLMax; liSCL++, riSCL++, tiTCL++)
        {
          val += *riSCL - *liSCL;
          (*tiTCL) = val*iarr;
        }

        for(glm::vec3* tiTCLMax=tiTCL+nTail; tiTCL<tiTCLMax; liSCL++, tiTCL++)
        {
          val += lv - *liSCL;
          (*tiTCL) = val*iarr;
        }
      }
    }
  });
}

//##################################################################################################
void boxBlurT_4(std::vector<glm::vec3>& scl, std::vector<glm::vec3>& tcl, size_t w, size_t h, size_t r)
{
  float iarr = 1.0f / float(r + r + 1);

  size_t c=0;
  tp_utils::parallel([&](auto locker)
  {
    for(;;)
    {
      size_t i;
      locker([&]{i=c; c++;});

      if(i>=w)
        return;

      size_t ti = i;
      size_t li = ti;
      size_t ri = ti + r * w;

      glm::vec3 fv = scl[ti];
      glm::vec3 lv = scl[ti + w * (h - 1)];
      glm::vec3 val = float(r + 1)*fv;

      for(size_t j=0; j<r; j++)
        val += scl[ti + j * w];

      for(size_t j=0; j<=r; j++)
      {
        val += scl[ri] - fv;
        tcl[ti] = val*iarr;
        ri += w; ti += w;
      }

      for(size_t j=r+1; j<h-r; j++)
      {
        val += scl[ri] - scl[li];
        tcl[ti] = val*iarr;
        li += w;
        ri += w;
        ti += w;
      }

      for(size_t j=h-r; j<h; j++)
      {
        val += lv-scl[li];
        tcl[ti] = val*iarr;
        li += w;
        ti += w;
      }
    }
  });
}

//##################################################################################################
void boxBlur_4(std::vector<glm::vec3>& scl, std::vector<glm::vec3>& tcl, size_t w, size_t h, size_t r)
{
  std::memcpy(tcl.data(), scl.data(), sizeof(glm::vec3)*scl.size());
  boxBlurH_4(tcl, scl, w, h, r);
  boxBlurT_4(scl, tcl, w, h, r);
}

//##################################################################################################
void gaussBlur_4(std::vector<glm::vec3>& scl, std::vector<glm::vec3>& tcl, size_t w, size_t h, size_t r)
{
  std::vector<size_t> bxs = boxesForGauss(float(r), 3);
  boxBlur_4(scl, tcl, w, h, (bxs[0] - 1) / 2);
  boxBlur_4(tcl, scl, w, h, (bxs[1] - 1) / 2);
  boxBlur_4(scl, tcl, w, h, (bxs[2] - 1) / 2);
}

}



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
std::vector<float> ConvolutionMatrix::matrixDataF()const
{
  std::vector<float> r;
  r.reserve(m_matrixData.size());
  for(double v : m_matrixData)
    r.push_back(float(v));
  return r;
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

  divideBySize();
}

//##################################################################################################
void ConvolutionMatrix::makeBlur3()
{
  m_width =  3;
  m_height = 3;

  m_matrixData =
  {
    1, 3, 1,
    3, 5, 3,
    1, 3, 1
  };

  divideBySize();
}

//##################################################################################################
void ConvolutionMatrix::makeBlur5()
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

  divideBySize();
}

//##################################################################################################
void ConvolutionMatrix::divideBySize()
{
  double f = 1.0 / (double(m_width) * double(m_height));
  for(double& v : m_matrixData)
    v *= f;
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

  if((width*height)>matrixData.size())
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
        d->r = uint8_t(tpBound(0, int(s->red  ), 255));
        d->g = uint8_t(tpBound(0, int(s->green), 255));
        d->b = uint8_t(tpBound(0, int(s->blue ), 255));
        d->a = 255;
        d++;
        s++;
      }
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ColorMapF convolvePadded(const tp_image_utils::ColorMapF& src, const std::vector<float>& matrixData, size_t width, size_t height)
{
  if(width<=1 || height<=1)
    return tp_image_utils::ColorMapF();

  if((width*height)>matrixData.size())
    return tp_image_utils::ColorMapF();

  size_t dw = src.width()  - (width-1);
  size_t dh = src.height() - (height-1);

  size_t marginX = (width-1)/2;
  size_t marginY = (height-1)/2;

  if(dw<1 || dh<1 || dw>src.width() || dh>src.height())
    return tp_image_utils::ColorMapF();

  tp_image_utils::ColorMapF dst{src.width(), src.height(), nullptr, {0.0f,0.0f,0.0f,0.0f}};
  glm::vec4* bufferData = dst.data();
  for(size_t my=0; my<height; my++)
  {
    for(size_t mx=0; mx<width;  mx++)
    {
      float weight = matrixData.at((my*width)+mx);
      size_t dyCounter=0;
      tp_utils::parallel([&](auto locker)
      {
        for(;;)
        {
          size_t dy;
          locker([&]{dy=dyCounter; dyCounter++;});

          if(dy>=dh)
            return;

          auto* s = src.constData() + ((dy+my)*src.width()) + mx;
          glm::vec4* d = bufferData + ((dy+marginY)*dst.width()) + marginX;
          glm::vec4* dMax = d + dw;
          for(; d<dMax; d++, s++)
            (*d) += (*s) * weight;
        }
      });
    }
  }

  // Copy in the top and bottom margins
  for(size_t y1=0; y1<marginY; y1++)
  {
    size_t y2 = (dst.height()-1) - y1;
    std::memcpy(dst.data()+dst.width()*y1, src.constData()+src.width()*y1, src.width()*sizeof(glm::vec4));
    std::memcpy(dst.data()+dst.width()*y2, src.constData()+src.width()*y2, src.width()*sizeof(glm::vec4));
  }

  // Copy in the left and right
  {
    size_t xSize = marginX*sizeof(glm::vec4);
    size_t x1 = 0;
    size_t x2 = (dst.width()-1) - marginX;

    size_t yMax = dst.height() - marginY;

    for(size_t y=marginY; y<yMax; y++)
    {
      std::memcpy(dst.data()+(y*dst.width())+x1, src.constData()+(y*src.width())+x1, xSize);
      std::memcpy(dst.data()+(y*dst.width())+x2, src.constData()+(y*src.width())+x2, xSize);
    }
  }

  return dst;
}

//##################################################################################################
std::vector<glm::vec3> gaussBlur(std::vector<glm::vec3>& source,
                                 size_t w,
                                 size_t h,
                                 size_t radius)
{
  std::vector<glm::vec3> target(source.size(), glm::vec3(0.0f,0.0f,0.0f));
  gaussBlur_4(source, target, w, h, radius);
  return target;
}

//##################################################################################################
tp_image_utils::ColorMap gaussBlur(const tp_image_utils::ColorMap& source, size_t radius)
{
  std::vector<glm::vec3> in(source.size(), glm::vec3(0.0f,0.0f,0.0f));

  {
    auto d = in.data();
    auto s = source.constData();
    auto sMax = s+source.size();

    for(; s<sMax; s++, d++)
    {
      d->x = s->r;
      d->y = s->g;
      d->z = s->b;
    }
  }

  std::vector<glm::vec3> out = gaussBlur(in, source.width(), source.height(), radius);

  tp_image_utils::ColorMap result(source.width(), source.height());
  {
    auto d = result.data();
    auto s = out.data();
    auto sMax = s+out.size();

    for(; s<sMax; s++, d++)
    {
      d->r = uint8_t(s->x);
      d->g = uint8_t(s->y);
      d->b = uint8_t(s->z);
      d->a = 255;
    }
  }

  return result;
}

//##################################################################################################
tp_image_utils::ColorMapF gaussBlur(const tp_image_utils::ColorMapF& source, size_t radius)
{
  std::vector<glm::vec3> in(source.size(), glm::vec3(0.0f,0.0f,0.0f));

  {
    auto d = in.data();
    auto s = source.constData();
    auto sMax = s+source.size();

    for(; s<sMax; s++, d++)
    {
      d->x = s->x;
      d->y = s->y;
      d->z = s->z;
    }
  }

  std::vector<glm::vec3> out = gaussBlur(in, source.width(), source.height(), radius);

  tp_image_utils::ColorMapF result(source.width(), source.height());
  {
    auto d = result.data();
    auto s = out.data();
    auto sMax = s+out.size();

    for(; s<sMax; s++, d++)
    {
      d->x = s->x;
      d->y = s->y;
      d->z = s->z;
      d->w = 1.0f;
    }
  }

  return result;

}

//##################################################################################################
tp_image_utils::ColorMapF blur3(const tp_image_utils::ColorMapF& src)
{
  ConvolutionMatrix m;
  m.makeBlur3();
  return convolvePadded(src, m.matrixDataF(), m.width(), m.height());
}

//##################################################################################################
tp_image_utils::ColorMapF blur5(const tp_image_utils::ColorMapF& src)
{
  ConvolutionMatrix m;
  m.makeBlur5();
  return convolvePadded(src, m.matrixDataF(), m.width(), m.height());
}


}
