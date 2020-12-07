#include "tp_image_utils_functions/SlotFill.h"
#include "tp_image_utils_functions/SignedDistanceField.h"

#include "tp_utils/TimeUtils.h"
#include "tp_utils/DebugUtils.h"

#include <cmath>

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ByteMap slotFill(const tp_image_utils::ByteMap& src, const SlotFillParameters& params)
{
  const uint8_t solid = params.solid;
  const uint8_t slot  = params.slot;

  tp_image_utils::ByteMap dst = src;

  float hw = float(dst.width())/1.0f;
  float hh = float(dst.height())/1.0f;
  float fw = float(dst.width());
  float fh = float(dst.height());

  size_t maxAngle = params.startAngle + params.maxAngle;
  for(size_t a=params.startAngle; a<=maxAngle; a+=params.stepAngle)
  {
    float af = float(a) * static_cast<float>(0.01745329251994329576923690768489);
    float saf = std::sin(af);
    float caf = std::cos(af);

    //Calculate how large the rotated image will be
    tp_image_utils::ByteMap rot;
    {
      float ix = std::fabs(saf*fw + caf*fh);
      float iy = std::fabs(caf*fw + saf*fh);
      float jx = std::fabs(saf*fw);
      float jy = std::fabs(caf*fw);

      float i = tpMax(ix, jx);
      float j = tpMax(iy, jy);

      rot = tp_image_utils::ByteMap(size_t(i), size_t(j));
    }
    size_t w = rot.width();
    size_t h = rot.height();

    float shw = float(rot.width())/1.0f;
    float shh = float(rot.height())/1.0f;

    auto proj = [=](float sx, float sy, size_t& dx, size_t& dy)
    {
      sx -= shw;
      sy -= shh;

      dx = size_t((saf*sx + caf*sy) + hw);
      dy = size_t((caf*sx + saf*sy) + hh);
    };

    //Rotate the image
    {
      uint8_t* r = rot.data();
      for(size_t y=0; y<h; y++)
      {
        for(size_t x=0; x<w; x++, r++)
        {
          size_t px=0;
          size_t py=0;
          proj(float(x), float(y), px, py);
          (*r) = (px<src.width() && py<src.height())?src.pixel(px, py) : params.slot;
        }
      }
    }


    const uint8_t* rm1 = rot.constData();
    const uint8_t* r = rm1+w;
    const uint8_t* rp1 = r+w;
    size_t hm1 = h-1;
    for(size_t y=1; y<hm1; y++)
    {
      size_t startFrame = 0;
      size_t slotFrame  = 0;
      size_t endFrame   = 0;

      size_t startFrameBorderM1=0;
      size_t startFrameBorderP1=0;

      size_t x=0;
      auto incAll = [&]()
      {
        x++;
        rm1++;
        r++;
        rp1++;
      };

      //Find the first pixel of the first startFrame
      for(; x<w; incAll())
      {
        if((*r) == solid)
        {
          startFrame = 1;
          incAll();
          break;
        }
      }

      //Find the last pixel of a startFrame
      size_t slotStart  = 0; //The first pixel of a slot
      for(; x<rot.width(); incAll())
      {
        if((*r) == slot)
        {
          slotStart = x;
          incAll();
          break;
        }
        startFrame++;
        startFrameBorderM1 = ((*rm1)==slot)?startFrameBorderM1+1:0;
        startFrameBorderP1 = ((*rp1)==slot)?startFrameBorderP1+1:0;
      }

      while(x<rot.width())
      {
        size_t slotEnd    = 0; //The last pixel of a slot

        //Find the last pixel of a slotFrame
        for(; x<rot.width(); incAll())
        {
          if((*r) == solid)
          {
            slotEnd = x;
            incAll();
            break;
          }
          slotFrame++;
        }

        //Find the last pixel of a endFrame
        size_t nextSlotStart = 0;
        size_t endFrameBorderM1=0; bool endFrameBorderM1Done=false;
        size_t endFrameBorderP1=0; bool endFrameBorderP1Done=false;
        size_t startFrameBorderM1Next=0;
        size_t startFrameBorderP1Next=0;
        for(; x<rot.width(); incAll())
        {
          if((*r) == slot)
          {
            nextSlotStart = x;
            incAll();
            break;
          }
          endFrame++;

          if(!endFrameBorderM1Done)
          {
            if((*rm1)==slot)
              endFrameBorderM1++;
            else
              endFrameBorderM1Done=true;
          }

          if(!endFrameBorderP1Done)
          {
            if((*rp1)==slot)
              endFrameBorderP1++;
            else
              endFrameBorderP1Done=true;
          }

          startFrameBorderM1Next = ((*rm1)==slot)?startFrameBorderM1Next+1:0;
          startFrameBorderP1Next = ((*rp1)==slot)?startFrameBorderP1Next+1:0;
        }

        size_t startBorder = tpMax(startFrameBorderM1, startFrameBorderP1);
        size_t   endBorder = tpMax(  endFrameBorderM1,   endFrameBorderP1);

        if(slotFrame>0                                  &&

           startFrame>=params.endMinEach                &&
           endFrame>=params.endMinEach                  &&
           (startFrame>=params.endMinEither          ||
            endFrame>=params.endMinEither)              &&
           (startFrame+endFrame)>=params.endMinSum      &&

           startBorder>=params.borderMinEach            &&
           endBorder>=params.borderMinEach              &&
           (startBorder>=params.borderMinEither      ||
            endBorder>=params.borderMinEither)          &&
           (startBorder+endBorder)>=params.borderMinSum)
        {
          for(size_t i=slotStart; i<=slotEnd; i++)
          {
            size_t px=0;
            size_t py=0;
            proj(float(i), float(y), px, py);

            if(px<src.width() && py<src.height())
              dst.setPixel(px, py, solid);
          }
        }

        startFrame = endFrame;
        slotFrame  = 0;
        endFrame   = 0;
        slotStart  = nextSlotStart;
        startFrameBorderM1 = startFrameBorderM1Next;
        startFrameBorderP1 = startFrameBorderP1Next;
      }
    }
  }

  return dst;
}

}
