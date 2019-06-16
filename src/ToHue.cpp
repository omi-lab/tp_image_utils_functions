#include "tp_image_utils_functions/ToHue.h"

#include "tp_image_utils/ColorMap.h"

#include <cmath>

namespace tp_image_utils_functions
{
//##################################################################################################
tp_image_utils::ColorMap toHue(const tp_image_utils::ColorMap& src)
{
  tp_image_utils::ColorMap dst(src.width(), src.height());

  TPPixel* d = dst.data();
  const TPPixel* s = src.constData();
  const TPPixel* sMax = s + src.size();

  while(s<sMax)
  {
    int r = s->r;
    int g = s->g;
    int b = s->b;

    int t = r + g + b;

    if(t==0 || t==765)
    {
      r = 85;
      g = 85;
      b = 85;
    }
    else
    {
      float f = 255.0f / float(t);
      r = int(float(r) * f);
      g = int(float(g) * f);
      b = int(float(b) * f);
    }

    d->r = uint8_t(r);
    d->g = uint8_t(g);
    d->b = uint8_t(b);
    d->a = 255;

    s++;
    d++;
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ColorMap toHue(const tp_image_utils::ByteMap& src)
{
  tp_image_utils::ColorMap dst(src.width(), src.height());

  const uint8_t* s = src.constData();
  const uint8_t* sMax = s + src.size();
  TPPixel* d = dst.data();

  while(s<sMax)
  {
    float h = float(*s) * 1.411f;
    h = h * 180.0f / 3.14f;

    d->r = uint8_t(tpBound(0, int(std::sin(h) * 255.0f), 255));
    d->g = 0;//uint8_t((1.0f+sin(h+2.0943951f)) * 1.0f);
    d->b = 0;//uint8_t((1.0f+sin(h+4.1887902f)) * 1.0f);
    d->a = 255;

    s++;
    d++;
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap toHueGray(const tp_image_utils::ColorMap& src)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());

  uint8_t* d = dst.data();
  const TPPixel* s = src.constData();
  const TPPixel* sMax = s + src.size();

    while(s<sMax)
    {
      float r = s->r;
      float g = s->g;
      float b = s->b;

      float min = tpMin(tpMin(r, g), b);
      float max = tpMax(tpMax(r, g), b);

      float hue = 0.0f;
      if (std::fabs(max-r)<1.0f)
        hue = (g - b) / (max - min);

      else if (std::fabs(max-g)<1.0f)
        hue = 2.0f + (b - r) / (max - min);

      else
        hue = 4.0f + (r - g) / (max - min);


      hue = hue * 42.5f;
      if (hue < 0) hue = hue + 255.0f;

      float val = (hue>0.0f)?(hue+0.5f):(hue-0.5f);

      (*d) = uint8_t(val);

      s++;
      d++;
    }

  return dst;
}
}
