#include "tp_image_utils_functions/NormalizeBrightness.h"

#include "tp_utils/DebugUtils.h"

#include <array>

namespace tp_image_utils_functions
{

//##################################################################################################
const char* shiftBrightnessModeToString(ShiftBrightnessMode mode)
{
  switch(mode)
  {
    case ShiftBrightnessMode::None       : return "None";
    case ShiftBrightnessMode::ByValue    : return "By value";
    case ShiftBrightnessMode::ByMean     : return "By mean";
    case ShiftBrightnessMode::ByMode     : return "By mode";
    case ShiftBrightnessMode::ByMedian   : return "By median";
    case ShiftBrightnessMode::BySoftMode : return "By soft mode";
  }
  return "None";
}

//##################################################################################################
ShiftBrightnessMode shiftBrightnessModeFromString(const std::string& mode)
{
  if(mode == "By value")     return ShiftBrightnessMode::ByValue;
  if(mode == "By mean")      return ShiftBrightnessMode::ByMean;
  if(mode == "By mode")      return ShiftBrightnessMode::ByMode;
  if(mode == "By median")    return ShiftBrightnessMode::ByMedian;
  if(mode == "By soft mode") return ShiftBrightnessMode::BySoftMode;

  return ShiftBrightnessMode::None;
}

//##################################################################################################
void normalizeBrightness(tp_image_utils::ColorMap& image, int radius, NormalizationMode mode, float exaggeration)
{
  if(mode == NormalizationMode::None)
    return;

  float exaggerationI = 1.0f / exaggeration;

  int w = image.width();
  int h = image.height();

  tp_image_utils::ColorMap results(w, h);

  int stride = image.width();

  const TPPixel* src = image.constData();
  TPPixel* dst = results.data();

  int xMax = (image.width());
  int yMax = image.height();
  for(int y=0; y<yMax; y++)
  {
    int ys = (y*stride);

    int ryMin = tpMax(0, (y-radius));
    int ryMax = tpMin(yMax, (y+radius));
    for(int x=0; x<xMax; x++)
    {
      int brightness=0;
      int count=0;
      for(int py=ryMin; py<ryMax; py++)
      {
        int pys = (py*stride);
        int ry = abs(y-py)/2;

        int rxMin = tpMax(0, (x-((radius-ry))));
        int rxMax = tpMin(xMax, ((x+(radius-ry))));
        for(int px=rxMin; px<rxMax; px++)
        {
          const TPPixel& p = src[pys+px];
          brightness += int(p.r);
          brightness += int(p.g);
          brightness += int(p.b);
          count+=3;
        }
      }

      const TPPixel& ps = src[ys+x];
      float bluef  = float(ps.b)/255.0f;
      float greenf = float(ps.g)/255.0f;
      float redf   = float(ps.r)/255.0f;

      float bf = 1.0f;
      if(mode == NormalizationMode::Normalize)
        bf     = (1.5f-((float(brightness)/float(count))/255.0f));
      else
      {
        float av  = (float(brightness)/float(count))/255.0f;
        float pav = (bluef+greenf+redf)/3.0f;
        bf        = (pav>av)?exaggeration:exaggerationI;
      }

      TPPixel& pd = dst[ys+x];
      pd.b = uint8_t(tpMin(255.0f, 255.0f*( bluef*bf)));
      pd.g = uint8_t(tpMin(255.0f, 255.0f*(greenf*bf)));
      pd.r = uint8_t(tpMin(255.0f, 255.0f*(  redf*bf)));
      pd.a = 255;
    }
  }

  image = results;
}

//##################################################################################################
void shiftBrightness(tp_image_utils::ColorMap& image, ShiftBrightnessMode mode, uint8_t value)
{
  if(image.size()<1)
    return;

  size_t shift=0;
  {
    const TPPixel* s = image.constData();
    const TPPixel* sMax = s+image.size();
    switch(mode)
    {
    case ShiftBrightnessMode::None:
    {
      return;
    }

    case ShiftBrightnessMode::ByValue:
    {
      shift = value;
      break;
    }

    case ShiftBrightnessMode::ByMean:
    {
      uint64_t total=0;
      uint64_t count=image.size()*3;
      for(; s<sMax; s++)
        total+=s->r+s->g+s->b;
      shift = total / count;
      break;
    }

    case ShiftBrightnessMode::ByMode:
    {
      std::array<int, 256> bins{0};
      //bins.fill(0);
      for(; s<sMax; s++)
        bins.at(size_t(s->r+s->g+s->b)/3)++;
      int max = 0;
      for(size_t i=0; i<bins.size(); i++)
      {
        if(bins.at(i)>max)
        {
          max = bins.at(i);
          shift = i;
        }
      }
      break;
    }

    case ShiftBrightnessMode::ByMedian:
    {
      std::vector<int> list;
      list.reserve(image.size());
      for(; s<sMax; s++)
        list.push_back(s->r+s->g+s->b);
      std::sort(list.begin(), list.end());
      shift = list[list.size()/2]/3;
      break;
    }

    case ShiftBrightnessMode::BySoftMode:
    {
      std::array<int, 26> bins{0};
      //bins.fill(0);
      for(; s<sMax; s++)
        bins.at(int(s->r+s->g+s->b)/30)++;
      int max = 0;
      for(size_t i=0; i<26; i++)
      {
        if(bins.at(i)>max)
        {
          max = bins.at(i);
          shift = i*10;
        }
      }
      break;
    }
    }

    if(mode != ShiftBrightnessMode::ByValue)
      shift = 128 - shift;
  }

  {
    TPPixel* s = image.data();
    TPPixel* sMax = s+image.size();
    for(; s<sMax; s++)
    {
      s->r = tpBound(0, int(s->r+shift), 255);
      s->g = tpBound(0, int(s->g+shift), 255);
      s->b = tpBound(0, int(s->b+shift), 255);
    }
  }
}

}
