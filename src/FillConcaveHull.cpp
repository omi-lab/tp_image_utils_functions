#include "tp_image_utils_functions/FillConcaveHull.h"
#include "tp_image_utils_functions/SignedDistanceField.h"

#include "tp_utils/TimeUtils.h"
#include "tp_utils/DebugUtils.h"

#include <cmath>

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ByteMap fillConcaveHull(const tp_image_utils::ByteMap& src, const FillConcaveHullParameters& params)
{
  tp_utils::ElapsedTimer t;
  t.start();

  size_t w = src.width();
  size_t h = src.height();

  if(w<3 || h<3)
    return src;

  struct Pixel
  {
    int32_t value;
    uint32_t direction;
  };

  //-- Create and populate the scratch buffer ------------------------------------------------------
  size_t firstH=w-1;
  size_t firstV=h-1;
  size_t lastH=1;
  size_t lastV=1;

  std::vector<Pixel> buffer;
  buffer.resize(src.size());
  {
    const uint8_t* s = src.constData();
    const uint8_t* sMax = s + src.size();
    Pixel* d = buffer.data();
    for(size_t y=0; s<sMax; y++)
    {
      for(size_t x=0; x<w; s++, d++, x++)
      {
        if((*s) == params.solid)
        {
          if(x<firstH)
            firstH=x;
          if(y<firstV)
            firstV=y;
          if(x>lastH)
            lastH=x;
          if(y>lastV)
            lastV=y;
          d->value = 0;
        }
        else
          d->value = 1;

        d->direction = 0;
      }
    }
  }

  auto getVal=[&](size_t x, size_t y)
  {
    return buffer.at((y*w)+x).value;
  };

  const auto& getPixel=[&](size_t x, size_t y)
  {
    return buffer.at((y*w)+x);
  };

  auto setVal=[&](size_t x, size_t y, int v, uint32_t d)
  {
    auto& b = buffer[(y*w)+x];
    b.value = v;
    b.direction = d;
  };

  size_t sxMax = tpMin(lastH+1, w-2);
  size_t syMax = tpMin(lastV+1, h-2);
  size_t sxMin = tpMax(size_t(2), firstH-1);
  size_t syMin = tpMax(size_t(2), firstV-1);

  int label=1;

  struct Vec
  {
    int x;
    int y;
    Vec(int x_, int y_):x(x_),y(y_){}
  };
  const std::array<Vec, 4> vectorsA = {{{0,-1},{-1,0},{0,1},{1,0}}};
  const std::array<Vec, 4> vectorsB = {{{1,0},{0,1},{-1,0},{0,-1}}};
  auto vectors = &vectorsA;

  uint32_t currentVector=0;
  auto advance = [&](size_t& x, size_t& y)
  {
    for(size_t i=0; i<2; i++)
    {
      Vec v = vectors->at(currentVector);

      if(getVal(x+size_t(v.x), y+size_t(v.y))!=0)
      {
        x+=size_t(v.x);
        y+=size_t(v.y);
        return true;
      }

      currentVector = (currentVector+1)&3;
    }

    return false;
  };

  bool workDone = false;
  int lastWorkLabel=2;
  do
  {
    workDone = false;
    for(size_t sy=syMin; sy<syMax; sy++)
    {
      for(size_t sx=sxMin; sx<sxMax; sx++)
      {
        if(getVal(sx, sy) == 0)
          continue;

        if(getVal(sx, sy-1)!=0 &&
           getVal(sx-1, sy)!=0 &&
           getVal(sx+1, sy)!=0 &&
           getVal(sx, sy+1)!=0)
          continue;

        bool workDoneInPass=false;
        for(int passc=0; passc<2 && !workDoneInPass; passc++)
        {
          vectors = (vectors == &vectorsA)?&vectorsB:&vectorsA;

          for(size_t i=0; i<4; i++)
          {
            Vec v1 = vectors->at(currentVector);
            currentVector = (currentVector+1)&3;
            Vec v2 = vectors->at(currentVector);

            if(getVal(sx+size_t(v1.x), sy+size_t(v1.y))==0 &&
               getVal(sx+size_t(v2.x), sy+size_t(v2.y))!=0)
              break;
          }

          label++;

          size_t cx = sx;
          size_t cy = sy;

          //-- Trace with label --------------------------------------------------------------------
          setVal(cx, cy, label, currentVector);
          bool continueToFill = true;
          for(;;)
          {
            if(!advance(cx, cy))
            {
              setVal(cx, cy, 0, 0);
              lastWorkLabel = label+1;
              workDone = true;
              workDoneInPass = true;
              continueToFill = false;
              break;
            }

            const Pixel& cp = getPixel(cx, cy);
            if(cp.value == label)
            {
              currentVector = cp.direction;
              lastWorkLabel = label+1;
              workDone = true;
              workDoneInPass = true;
              break;
            }

            if(cp.value>lastWorkLabel && cp.direction==currentVector)
            {
              continueToFill = false;
              break;
            }

            if(cx<=sxMin || cy<=syMin || cx>=sxMax || cy>=syMax)
            {
              continueToFill = false;
              break;
            }

            if((getVal(cx, cy-1)==0 && getVal(cx, cy+1)==0) ||
               (getVal(cx-1, cy)==0 && getVal(cx+1, cy)==0))
            {
              lastWorkLabel = label+1;
              workDone = true;
              workDoneInPass = true;
              setVal(cx, cy, 0, currentVector);
            }
            else
              setVal(cx, cy, label, currentVector);
          }

          if(continueToFill)
          {
            do
            {
              setVal(cx, cy, 0, 0);
            }while(advance(cx, cy));
          }
        }
      }
    }
  }while(workDone);

  tp_image_utils::ByteMap dst(w, h);
  {
    uint8_t slot = 255 - params.solid;
    const Pixel* s = buffer.data();
    const Pixel* sMax = s + buffer.size();
    uint8_t* d = dst.data();
    for(; s<sMax; s++, d++)
      (*d) = (s->value==0)?params.solid:slot;
  }

  return dst;
}

}
