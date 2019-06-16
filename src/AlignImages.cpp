#include "tp_image_utils_functions/AlignImages.h"
#include "tp_image_utils_functions/ExtractRect.h"

#include <cmath>

namespace tp_image_utils_functions
{

//##################################################################################################
AlignImages::Rect AlignImages::Rect::intersected(const AlignImages::Rect& other)
{
  int l = tpMax(x,   other.x);
  int r = tpMin(x+w, other.x+other.w);

  int t = tpMax(y,   other.y);
  int b = tpMin(y+h, other.y+other.h);

  if(l>=r)
    return Rect();

  if(t>=b)
    return Rect();

  return Rect(l, t, r-l, b-t);
}

//##################################################################################################
std::pair<size_t, size_t> AlignImages::calculateMicroAlignment(const tp_image_utils::ByteMap& reference,
                                                               const tp_image_utils::ByteMap& other,
                                                               size_t maxOffset)
{
  std::pair<size_t, size_t> translation(0, 0);

  size_t steps = 100;
  size_t margin = maxOffset+3;
  size_t minD = steps + (margin*2) + 2;
  size_t best = (steps*steps) * 255;


  if(reference.size() != other.size())
    return std::pair<int, int>(0, 0);

  if(reference.width()>=minD && reference.height()>=minD)
  {
    size_t sx = (reference.width()-(2*margin)) / steps;
    size_t sy = (reference.height()-(2*margin)) / steps;

    for(size_t x=-maxOffset; x<=maxOffset; x++)
    {
      for(size_t y=-maxOffset; y<=maxOffset; y++)
      {
        size_t s = 0;
        for(size_t i=0; i<steps; i++)
        {
          for(size_t j=0; j<steps; j++)
          {
            size_t px = (sx * i) + margin;
            size_t py = (sy * j) + margin;
            s += size_t(std::abs(int(reference.pixel(px, py)) - int(other.pixel(px-x, py-y))));
          }
        }

        if(s<best)
        {
          best=s;
          translation.first  = x;
          translation.second = y;
        }
      }
    }
  }
  return translation;
}

//##################################################################################################
void AlignImages::translateAndClipPair(const std::pair<int, int>& translation,
                                       tp_image_utils::ColorMap& reference,
                                       tp_image_utils::ColorMap& image)
{
  int tx = translation.first;
  int ty = translation.second;

  Rect refRect(tx, ty, int(reference.width()), int(reference.height()));
  Rect imgRect(tx, ty, int(image.width()), int(image.height()));

  Rect i = refRect.intersected(imgRect);

  reference = reference.subImage(size_t(i.x), size_t(i.y), size_t(i.w)+size_t(i.w), size_t(i.y)+size_t(i.h));
  i.x-=tx;
  i.y-=ty;
  image = image.subImage(size_t(i.x), size_t(i.y), size_t(i.w)+size_t(i.w), size_t(i.y)+size_t(i.h));
}

//##################################################################################################
AlignImages::SkewedRegion AlignImages::calculateSkewedRegion(const tp_image_utils::ByteMap& reference,
                                                             const tp_image_utils::ByteMap& other,
                                                             size_t maxOffset)
{
  AlignImages::SkewedRegion result;

  result.referenceRect.w = int(reference.width());
  result.referenceRect.h = int(reference.height());

  const int steps = 100;
  size_t margin = maxOffset+3;
  size_t minD = steps + (margin*2) + 2;

  if(reference.size() != other.size() || reference.width()<minD || reference.height()<minD)
    return result;

  struct CornerDetails_lt
  {
    std::pair<int, int> translation{std::pair<int, int>(0, 0)};
    float best{(steps*steps) * 255};
    float fx;
    float fy;

    float s{0.0f};
    float ax{0.0f};
    float ay{0.0f};
    float wt{0.0f};
    CornerDetails_lt(float fx_, float fy_):fx(fx_),fy(fy_){}
  };

  CornerDetails_lt topLeft    (0.0f, 0.0f);
  CornerDetails_lt topRight   (1.0f, 0.0f);
  CornerDetails_lt bottomRight(1.0f, 1.0f);
  CornerDetails_lt bottomLeft (0.0f, 1.0f);

  CornerDetails_lt* corners[]={&topLeft, &topRight, &bottomRight, &bottomLeft};

  size_t sx = (reference.width()-(2*margin)) / steps;
  size_t sy = (reference.height()-(2*margin)) / steps;

  float fw = reference.width();
  float fh = reference.height();

  for(int x=-int(maxOffset); x<=int(maxOffset); x++)
  {
    for(int y=-int(maxOffset); y<=int(maxOffset); y++)
    {
      for(CornerDetails_lt* cc : corners)
      {
        cc->s =0.0f;
        cc->ax=0.0f;
        cc->ay=0.0f;
        cc->wt=0.0f;
      }

      for(size_t i=0; i<steps; i++)
      {
        for(size_t j=0; j<steps; j++)
        {
          size_t px = (sx * i) + margin;
          size_t py = (sy * j) + margin;

          int opx = int(px)+x;
          int opy = int(py)+y;

          float ofx = float(opx) / fw;
          float ofy = float(opy) / fh;

          float s = abs(int(reference.pixel(px, py)) - int(other.pixel(size_t(opx), size_t(opy))));

          for(CornerDetails_lt* cc : corners)
          {
            float wt = 1.0f-(std::fabs(cc->fx-ofx) * std::fabs(cc->fy-ofy));

            cc->s  += wt*s;
            cc->ax += wt*float(x);
            cc->ay += wt*float(y);
            cc->wt += wt;
          }
        }
      }

      for(CornerDetails_lt* cc : corners)
      {
        float sw = cc->s/cc->wt;
        if(sw<cc->best && cc->wt>1.0f)
        {
          cc->best=sw;
          cc->translation.first  = int(cc->ax/cc->wt);
          cc->translation.second = int(cc->ay/cc->wt);
        }
      }
    }
  }

  for(CornerDetails_lt* cc : corners)
  {
    int px = int((cc->fx*fw) + float(cc->translation.first));
    int py = int((cc->fy*fh) + float(cc->translation.second));
    result.otherRegion.emplace_back(tp_image_utils::PointTypeRectCorner, px, py);
  }

  return result;
}

//##################################################################################################
void AlignImages::extractAndClipPair(const AlignImages::SkewedRegion& skewedRegion,
                                     tp_image_utils::ColorMap& reference,
                                     tp_image_utils::ColorMap& image)
{
  if(skewedRegion.referenceRect.x != 0                 ||
     skewedRegion.referenceRect.y != 0                 ||
     skewedRegion.referenceRect.w != int(reference.width()) ||
     skewedRegion.referenceRect.h != int(reference.height()))
  {
    const auto& r = skewedRegion.referenceRect;
    reference = reference.subImage(size_t(r.x),
                                   size_t(r.y),
                                   size_t(r.x)+size_t(r.w),
                                   size_t(r.y)+size_t(r.h));
  }

  if((skewedRegion.otherRegion.size() == 4)                                                 &&
     (skewedRegion.otherRegion.at(0).positionEquals(tp_image_utils::Point(0.0f,          0.0f           ))) &&
     (skewedRegion.otherRegion.at(1).positionEquals(tp_image_utils::Point(image.width(), 0.0f           ))) &&
     (skewedRegion.otherRegion.at(2).positionEquals(tp_image_utils::Point(image.width(), image.height() ))) &&
     (skewedRegion.otherRegion.at(3).positionEquals(tp_image_utils::Point(0.0f,          image.height() ))))
    return;

  std::vector<std::string> errors;
  image = ExtractRect::extractRect(image,
                                   skewedRegion.otherRegion,
                                   size_t(skewedRegion.referenceRect.w),
                                   size_t(skewedRegion.referenceRect.h),
                                   errors);
}

}
