#include "tp_image_utils_functions/ExtractRect.h"

namespace tp_image_utils_functions
{

namespace
{
//##################################################################################################
tp_image_utils::Point calculatePoint(const std::vector<tp_image_utils::Point>& points,
                                     float f)
{
  if(points.empty())
    return tp_image_utils::Point();

  size_t m = points.size()-1;

  f = f*float(m);

  size_t a = tpBound(size_t(0), size_t(f), m);
  size_t b = tpBound(size_t(0), a+1, m);

  if(a==b)
    return points.at(a);

  float af = float(b) - f;
  float bf = f - float(a);

  const auto& ap = points.at(a);
  const auto& bp = points.at(b);

  return tp_image_utils::Point(tp_image_utils::PointTypeGeneral, (ap.x*af)+(bp.x*bf), (ap.y*af)+(bp.y*bf));
}

//##################################################################################################
tp_image_utils::Point calculatePoint(const std::vector<tp_image_utils::Point>& topPoints,
                                     const std::vector<tp_image_utils::Point>& rightPoints,
                                     const std::vector<tp_image_utils::Point>& bottomPoints,
                                     const std::vector<tp_image_utils::Point>& leftPoints,
                                     float fx,
                                     float fy)
{
  tp_image_utils::Point top    = calculatePoint(topPoints,    fx);
  tp_image_utils::Point bottom = calculatePoint(bottomPoints, fx);
  tp_image_utils::Point left   = calculatePoint(leftPoints,   fy);
  tp_image_utils::Point right  = calculatePoint(rightPoints,  fy);

  return tp_image_utils::Point(tp_image_utils::PointTypeGeneral,
                               ((top   .x*(1.0f-fy))+
                                (bottom.x*      fy )+
                                (left  .x*(1.0f-fx))+
                                (right .x*      fx ))/2.0f,
                               ((top   .y*(1.0f-fy))+
                                (bottom.y*fy       )+
                                (left  .y*(1.0f-fx))+
                                (right .y*fx       ))/2.0f);
}
}

//##################################################################################################
tp_image_utils::ColorMap ExtractRect::extractRect(const tp_image_utils::ColorMap& sourceImage,
                                                  const std::vector<tp_image_utils::Point>& sourcePoints,
                                                  size_t width,
                                                  size_t height,
                                                  std::vector<std::string>& errors)
{
  size_t sw = sourceImage.width();
  size_t sh = sourceImage.height();

  if(width<1 || height<1 || sw<1 || sh<1)
    return tp_image_utils::ColorMap();

  sw--;
  sh--;

  std::vector<tp_image_utils::Point> topPoints;
  std::vector<tp_image_utils::Point> rightPoints;
  std::vector<tp_image_utils::Point> bottomPoints;
  std::vector<tp_image_utils::Point> leftPoints;

  if(sourcePoints.size()<4)
  {
    errors.emplace_back("To few source points!");
    return tp_image_utils::ColorMap();
  }

  topPoints.push_back(sourcePoints.at(0));
  size_t i=0;
  for(i++; i<sourcePoints.size(); i++)
  {
    const tp_image_utils::Point& point = sourcePoints.at(i);
    topPoints.push_back(point);

    if(i!=0 && point.type==tp_image_utils::PointTypeRectCorner)
    {
      rightPoints.push_back(point);
      break;
    }
  }

  for(i++; i<sourcePoints.size(); i++)
  {
    const tp_image_utils::Point& point = sourcePoints.at(i);
    rightPoints.push_back(point);

    if(point.type==tp_image_utils::PointTypeRectCorner)
    {
      bottomPoints.insert(bottomPoints.begin(), point);
      break;
    }
  }

  for(i++; i<sourcePoints.size(); i++)
  {
    const tp_image_utils::Point& point = sourcePoints.at(i);
    bottomPoints.insert(bottomPoints.begin(), point);

    if(point.type==tp_image_utils::PointTypeRectCorner)
    {
      leftPoints.insert(leftPoints.begin(), point);
      break;
    }
  }

  for(i++; i<sourcePoints.size(); i++)
  {
    const tp_image_utils::Point& point = sourcePoints.at(i);
    leftPoints.insert(leftPoints.begin(), point);

    if(point.type==tp_image_utils::PointTypeRectCorner)
    {
      errors.emplace_back("To many corners!");
      return tp_image_utils::ColorMap();
    }
  }

  if(leftPoints.empty())
  {
    errors.emplace_back("To few corners!");
    return tp_image_utils::ColorMap();
  }

  leftPoints.insert(leftPoints.begin(), sourcePoints.at(0));

  tp_image_utils::ColorMap result(width, height);

  for(size_t y=0; y<height; y++)
  {
    float fy = (float(y) / float(height));

    TPPixel* dst = result.data() + (y*width);
    for(size_t x=0; x<width; x++, dst++)
    {
      float fx = float(x) / float(width);

      tp_image_utils::Point ps = calculatePoint(topPoints, rightPoints, bottomPoints, leftPoints, fx, fy);

#if 0
      //Floor the coords and use that pixel
      size_t sx = tpBound(size_t(0), size_t(ps.x), sw);
      size_t sy = tpBound(size_t(0), size_t(ps.y), sh);

      (*dst) = sourceImage.pixel(sx, sy);
#endif

#if 1
      //Bilinear filtering
      size_t sx1 = tpBound(size_t(0), size_t(ps.x), sw);
      size_t sy1 = tpBound(size_t(0), size_t(ps.y), sh);
      size_t sx2 = tpMin(sx1+1, sw);
      size_t sy2 = tpMin(sy1+1, sh);

      float fxn = tpMin((ps.x-float(sx1)), 1.0f);
      float fyn = tpMin((ps.y-float(sy1)), 1.0f);
      float fxi = 1.0f-fxn;
      float fyi = 1.0f-fyn;

      TPPixel rgbTL = sourceImage.pixel(sx1, sy1);
      TPPixel rgbTR = sourceImage.pixel(sx2, sy1);
      TPPixel rgbBR = sourceImage.pixel(sx2, sy2);
      TPPixel rgbBL = sourceImage.pixel(sx1, sy2);

      float r = ((float(rgbTL.r)*fxi*fyi) + (float(rgbTR.r)*fxn*fyi) + (float(rgbBR.r)*fxn*fyn) + (float(rgbBL.r)*fxi*fyn));// / 2.0f;
      float g = ((float(rgbTL.g)*fxi*fyi) + (float(rgbTR.g)*fxn*fyi) + (float(rgbBR.g)*fxn*fyn) + (float(rgbBL.g)*fxi*fyn));// / 2.0f;
      float b = ((float(rgbTL.b)*fxi*fyi) + (float(rgbTR.b)*fxn*fyi) + (float(rgbBR.b)*fxn*fyn) + (float(rgbBL.b)*fxi*fyn));// / 2.0f;

      dst->r = uint8_t(r);
      dst->g = uint8_t(g);
      dst->b = uint8_t(b);
      dst->a = 255;
#endif

    }
  }

  return result;
}

//##################################################################################################
tp_image_utils::ColorMap ExtractRect::extractRect(const tp_image_utils::ColorMap& sourceImage,
                                                  const tp_image_utils::Grid& clippingGrid,
                                                  size_t width,
                                                  size_t height,
                                                  std::vector<std::string>& errors)
{
  std::vector<tp_image_utils::Point> sourcePoints;

  //Top left
  sourcePoints.push_back(clippingGrid.origin);

  //Top right
  {
    tp_image_utils::Point p = clippingGrid.origin;
    p.x += clippingGrid.xAxis.x*clippingGrid.xCells;
    p.y += clippingGrid.xAxis.y*clippingGrid.xCells;
    sourcePoints.push_back(p);
  }

  //Bottom right
  {
    tp_image_utils::Point p = clippingGrid.origin;
    p.x += clippingGrid.xAxis.x*clippingGrid.xCells;
    p.y += clippingGrid.xAxis.y*clippingGrid.xCells;
    p.x += clippingGrid.yAxis.x*clippingGrid.yCells;
    p.y += clippingGrid.yAxis.y*clippingGrid.yCells;
    sourcePoints.push_back(p);
  }

  //Bottom left
  {
    tp_image_utils::Point p = clippingGrid.origin;
    p.x += clippingGrid.yAxis.x*clippingGrid.yCells;
    p.y += clippingGrid.yAxis.y*clippingGrid.yCells;
    sourcePoints.push_back(p);
  }

  for(auto& i : sourcePoints)
    i.type = tp_image_utils::PointTypeRectCorner;

  return extractRect(sourceImage, sourcePoints, width, height, errors);
}

//##################################################################################################
tp_image_utils::ColorMap ExtractRect::extractRect(const tp_image_utils::ColorMap& sourceImage,
                                                  size_t x,
                                                  size_t y,
                                                  size_t w,
                                                  size_t h,
                                                  TPPixel pad)
{
  bool fill=false;

  size_t minSX=x;
  size_t maxSX=x+w;
  size_t minDX=0;
  if(maxSX>sourceImage.width())
  {
    maxSX=sourceImage.width();
    fill=true;
  }

  size_t minSY=y;
  size_t maxSY=y+h;
  size_t minDY=0;
  if(maxSY>sourceImage.height())
  {
    maxSY=sourceImage.height();
    fill=true;
  }

  tp_image_utils::ColorMap result(w, h);

  if(fill)
    result.fill(pad);

  size_t n = (maxSX - minSX) * sizeof(TPPixel);
  if(n>0)
  {
    size_t i=minDY;
    for(size_t j=minSY; j<maxSY; i++, j++)
    {
      const TPPixel* src = sourceImage.constData()+(j*sourceImage.width())+minSX;
      TPPixel* dst = result.data()+(i*result.width())+minDX;
      memcpy(dst, src, n);
    }
  }

  return result;
}

//##################################################################################################
tp_image_utils::ColorMap ExtractRect::extractRect(const tp_image_utils::ColorMap& sourceImage,
                                                  size_t x,
                                                  size_t y,
                                                  size_t w,
                                                  size_t h)
{
  size_t minSX=x;
  size_t maxSX=x+w;
  size_t minDX=0;
  if(maxSX>sourceImage.width())
    maxSX=sourceImage.width();

  size_t minSY=y;
  size_t maxSY=y+h;
  size_t minDY=0;
  if(maxSY>sourceImage.height())
    maxSY=sourceImage.height();

  tp_image_utils::ColorMap result(w, h);

  result.fill(TPPixel(0, 255, 0, 255));

  size_t copyW = (maxSX - minSX);
  size_t maxDX=minDX + copyW;
  size_t maxDY=minDY + (maxSY - minSY);
  size_t n = copyW * sizeof(TPPixel);
  if(n>0)
  {
    size_t i=0;
    size_t j=minSY;

    while(i<h)
    {
      const TPPixel* src = sourceImage.constData()+(j*sourceImage.width())+minSX;
      TPPixel* dst = result.data()+(i*result.width())+minDX;
      if(i<minDY || i>=maxDY)
      {
        const TPPixel* s = src;
        TPPixel* d = dst;
        TPPixel* dMax = d + copyW;
        for(; d<dMax; s++, d++)
        {
          (*d) = (*s);
          d->a = 1;
        }
      }
      else
        memcpy(dst, src, n);

      {
        TPPixel* dst2 = result.data();
        TPPixel pad = (*dst);
        pad.a = 1;
        for(size_t x=0; x<minDX; x++, dst2++)
          (*dst2) = pad;
      }

      {
        TPPixel* dst3 = result.data()+(i*result.width())+maxDX;
        TPPixel pad = *(dst3-1);
        pad.a = 1;
        for(size_t x=maxDX; x<w; x++, dst3++)
          (*dst3) = pad;
      }

      i++;
      if(i>minDY && j<(maxSY-1))
        j++;
    }
  }

  return result;
}

}
