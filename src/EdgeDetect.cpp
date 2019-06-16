#include "tp_image_utils_functions/EdgeDetect.h"

#include "tp_image_utils/ColorMap.h"

#include <cmath>

namespace tp_image_utils_functions
{

namespace
{

//##################################################################################################
int calcDist(int r1, int g1, int b1, int r2, int g2, int b2)
{
  int dr = r1-r2;
  int dg = g1-g2;
  int db = b1-b2;

  return int(std::sqrt(float(dr*dr + dg*dg + db*db)));
}

//##################################################################################################
int calcDist(int a, int b)
{
  int d = a-b;
  return abs(d);
}
}

//##################################################################################################
tp_image_utils::ByteMap edgeDetect(const tp_image_utils::ColorMap& src, int threshold)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());
  uint8_t* d = dst.data();

  int yMax = src.height()-1;
  for(int y=0; y<yMax; y++)
  {
    const TPPixel* s = src.constData() + (src.width()*y);
    const TPPixel* s2 = src.constData() + (src.width()*(y+1));
    const TPPixel* sMax = s + (src.width()-1);

    while(s<sMax)
    {
      if(calcDist(s->r, s->g, s->b, s[1].r, s[1].g, s[1].b)>threshold)
        *d = 255;
      else if(calcDist(s->r, s->g, s->b, s2[1].r, s2[1].g, s2[1].b)>threshold)
        *d = 255;
      else
        *d = 0;

      s++;
      s2++;
      d++;
    }

    *d = 0;
    d++;
  }

  int xMax = dst.width();
  for(int x=0; x<xMax; x++)
  {
    *d = 0;
    d++;
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap edgeDetect(const tp_image_utils::ByteMap& src, uint8_t threshold)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());
  uint8_t* d = dst.data();

  int yMax = src.height()-1;
  int xMax = src.width()-1;

  for(int y=0; y<yMax; y++)
  {
    for(int x=0; x<xMax; x++)
    {
      uint8_t c = src.pixel(x, y);
      if(calcDist(c, src.pixel(x+1, y))>threshold)
        *d = 255;
      else if(calcDist(c, src.pixel(x, y+1))>threshold)
        *d = 255;
      else
        *d = 0;

      d++;
    }

    *d = 0;
    d++;
  }

  xMax++;
  for(int x=0; x<xMax; x++)
  {
    *d = 0;
    d++;
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap edgeDetectCorner(const tp_image_utils::ByteMap& src, uint8_t threshold)
{
  int w = src.width();
  int h = src.height();

  tp_image_utils::ByteMap dst(w, h);
  uint8_t* d = dst.data();

  int yMax = src.height()-1;
  int xMax = src.width()-2;

  {
    uint8_t* dMax = d+w;
    for(; d<dMax; d++)
      *d = 0;
  }

  for(int y=1; y<yMax; y++)
  {
    *d = 0;

    const uint8_t* sm1 = src.constData() + (w*(y-1));
    const uint8_t* s   = sm1+w;
    const uint8_t* sp1 = s+w;

    d++;
    sm1++;
    s++;
    sp1++;

    uint8_t* dMax = d+xMax;
    for(; d<dMax; d++, sm1++, s++, sp1++)
    {
      uint8_t c = (*s);
      uint_fast8_t f =
          ((calcDist(c, *(s  -1))>threshold)?0:  1) +
          ((calcDist(c, *(sm1-1))>threshold)?0:  2) +
          ((calcDist(c, *(sm1  ))>threshold)?0:  4) +
          ((calcDist(c, *(sm1+1))>threshold)?0:  8) +
          ((calcDist(c, *(s  +1))>threshold)?0: 16) +
          ((calcDist(c, *(sp1+1))>threshold)?0: 32) +
          ((calcDist(c, *(sp1  ))>threshold)?0: 64) +
          ((calcDist(c, *(sp1-1))>threshold)?0:128);

      *d = ((f == 0b11100000) ||
            (f == 0b01110000) ||
            (f == 0b00111000) ||
            (f == 0b00011100) ||
            (f == 0b00001110) ||
            (f == 0b00000111) ||
            (f == 0b10000011) ||
            (f == 0b11000001) ||

            (f == 0b10100000) ||
            (f == 0b01010000) ||
            (f == 0b00101000) ||
            (f == 0b00010100) ||
            (f == 0b00001010) ||
            (f == 0b00000101) ||
            (f == 0b10000010) ||
            (f == 0b01000001))?255:0;
    }

    *d = 0;
    d++;
  }

  {
    uint8_t* dMax = d+w;
    for(; d<dMax; d++)
      *d = 0;
  }

  return dst;
}

}
