#include "tp_image_utils_functions/CellSegment.h"
#include "tp_image_utils_functions/SignedDistanceField.h"

#include "tp_image_utils/ColorMap.h"

#include "tp_utils/TimeUtils.h"
#include "tp_utils/DebugUtils.h"

#include <cmath>
#include <array>

namespace tp_image_utils_functions
{
namespace
{

//##################################################################################################
void floodGrowCell(tp_image_utils::ByteMap& result,
                   tp_image_utils::ByteMap& mask,
                   uint8_t cellID,
                   int x,
                   int y)
{
  size_t w = result.width();
  size_t h = result.height();

  tp_image_utils::ByteMap done(w, h);
  done.fill(0);
  std::vector<std::pair<size_t, size_t>> queue;
  queue.emplace_back(size_t(x), size_t(y));

  while(!queue.empty())
  {
    auto point = tpTakeLast(queue);
    size_t px = point.first;
    size_t py = point.second;

    if(px>=w || py>=h)
      continue;

    if(done.pixel(px, py))
      continue;
    done.setPixel(px, py, 1);

    if(mask.pixel(px, py) == 0)
      continue;

    mask.setPixel(px, py, 0);
    result.setPixel(px, py, cellID);

    queue.emplace_back(px+1, py  );
    queue.emplace_back(px,   py-1);
    queue.emplace_back(px-1, py  );
    queue.emplace_back(px  , py+1);
  }
}

//##################################################################################################
void boxGrowCell(tp_image_utils::ByteMap& result,
                 tp_image_utils::ByteMap& mask,
                 const CellSegmentParameters& params,
                 int w,
                 int h,
                 uint8_t cellID,
                 int v,
                 int x,
                 int y)
{
  int r    = (v*params.distanceFieldRadius)/256;// The radius of a circle that fits inside the square.
  int half = int(std::sqrt((r*r)/2));           // The half width of the square.

  size_t cxInt = size_t(tpMax(0, x-half));
  size_t cyInt = size_t(tpMax(0, y-half));
  size_t cxMax = size_t(tpMin(w, x+half));
  size_t cyMax = size_t(tpMin(h, y+half));

  bool canGrowXp = true; //Can grow up the X axis.
  bool canGrowYp = true; //Can grow up the Y axis.
  bool canGrowXm = true; //Can grow down the X axis.
  bool canGrowYm = true; //Can grow down the Y axis.

  do
  {
    do
    {
      if(canGrowXp)
      {
        bool& canGrow = canGrowXp;
        size_t& v     = cxMax;
        size_t  oInt  = cyInt;
        size_t  oMax  = cyMax;

        canGrow = false;

        size_t c = v+1;

        if(c>size_t(w))
          break;

        bool ok = true;
        for(size_t o=oInt; o<oMax; o++)
        {
          if(mask.pixel(c, o)<1)
          {
            ok = false;
            break;
          }
        }

        if(!ok)
          break;

        v = c;
        canGrow = true;
      }
    }
    while(false);

    do
    {
      if(canGrowYp)
      {
        bool&   canGrow = canGrowYp;
        size_t& v       = cyMax;
        size_t  oInt    = cxInt;
        size_t  oMax    = cxMax;

        canGrow = false;

        size_t c = size_t(v)+1;

        if(c>size_t(w))
          break;

        bool ok = true;
        for(size_t o=oInt; o<oMax; o++)
        {
          if(mask.pixel(o, c)<1)
          {
            ok = false;
            break;
          }
        }

        if(!ok)
          break;

        v = c;
        canGrow = true;
      }
    }
    while(false);

    do
    {
      if(canGrowXm)
      {
        bool&   canGrow = canGrowXm;
        size_t& v       = cxInt;
        size_t  oInt    = cyInt;
        size_t  oMax    = cyMax;

        canGrow = false;

        size_t c = v-1;

        if(c>=cxMax)
          break;

        bool ok = true;
        for(size_t o=oInt; o<oMax; o++)
        {
          if(mask.pixel(c, o)<1)
          {
            ok = false;
            break;
          }
        }

        if(!ok)
          break;

        v = c;
        canGrow = true;
      }
    }
    while(false);

    do
    {
      if(canGrowYm)
      {
        bool&   canGrow = canGrowYm;
        size_t& v       = cyInt;
        size_t  oInt    = cxInt;
        size_t  oMax    = cxMax;

        canGrow = false;

        size_t c = v-1;

        if(c>=cyMax)
          break;

        bool ok = true;
        for(size_t o=oInt; o<oMax; o++)
        {
          if(mask.pixel(o, c)<1)
          {
            ok = false;
            break;
          }
        }

        if(!ok)
          break;

        v = c;
        canGrow = true;
      }
    }
    while(false);
  }
  while(canGrowXp || canGrowYp || canGrowXm || canGrowYm);

  for(size_t cy=cyInt; cy<cyMax; cy++)
  {
    uint8_t* d    = mask.data() + (cy*size_t(w));
    uint8_t* dMax = d + cxMax;
    d += cxInt;
    for(; d<dMax; d++)
    {
      (*d) = 0;
    }
  }

  for(size_t cy=cyInt; cy<cyMax; cy++)
  {
    uint8_t* d    = result.data() + (cy*size_t(w));
    uint8_t* dMax = d + cxMax;
    d += cxInt;
    for(; d<dMax; d++)
    {
      (*d) = cellID;
    }
  }
}
}

//##################################################################################################
const char* cellGrowModeToString(CellGrowMode mode)
{
  switch(mode)
  {
  case CellGrowMode::Box:   return "Box";
  case CellGrowMode::Flood: return "Flood";
  }

  return "Box";
}

//##################################################################################################
CellGrowMode cellGrowModeFromString(const std::string& mode)
{
  if(mode == "Box")   return CellGrowMode::Box;
  if(mode == "Flood") return CellGrowMode::Flood;
  return CellGrowMode::Box;
}

//##################################################################################################
tp_image_utils::ByteMap cellSegmentInitialCells(const tp_image_utils::ByteMap& src,
                                                const CellSegmentParameters& params)
{
  size_t w = src.width();
  size_t h = src.height();

  tp_image_utils::ByteMap result(w, h);
  result.fill(0);

  //-- Create initial cells ------------------------------------------------------------------------
  //Here we try to fill the larger internal areas with squares, later we will pitch these squares
  //against each other to best fill the internal cavities.
  {
    //The mask will get filled out as we generate the initial cells
    tp_image_utils::ByteMap mask = src;
    uint8_t cellID = 0;
    for(int p=0; p<params.maxInitialCells; p++)
    {
      tp_utils::ElapsedTimer tt2;
      tt2.start();
      tp_image_utils::ByteMap ds = distanceField(mask, params.distanceFieldRadius);

      const uint8_t* s = ds.constData();
      const uint8_t* sMax = s + ds.size();
      uint8_t v=0;
      const uint8_t* b = nullptr;
      for(; s<sMax; s++)
      {
        if(v<(*s))
        {
          v = *s;
          b = s;
        }
      }

      if(b == nullptr)
        break;

      if(v<params.minRadius)
        break;

      cellID++;

      size_t i = size_t(b - ds.constData()); // Index of the most remote pixel.
      size_t y = i / w;                      // Y coordinate of the most remote pixel.
      size_t x = i - (y*w);                  // X coordinate of the most remote pixel.

      switch(params.cellGrowMode)
      {
      case CellGrowMode::Box:
        boxGrowCell(result, mask, params, int(w), int(h), cellID, v, int(x), int(y));
        break;

      case CellGrowMode::Flood:
        floodGrowCell(result, mask, cellID, int(x), int(y));
        break;
      }

      if(cellID==255)
        break;
    }
  }

  return result;
}

//##################################################################################################
tp_image_utils::ByteMap cellSegment(const tp_image_utils::ByteMap& src, const CellSegmentParameters& params)
{
  tp_image_utils::ByteMap labels = cellSegmentInitialCells(src, params);
  return cellSegment(src, labels, params);
}

//##################################################################################################
tp_image_utils::ByteMap cellSegment(const tp_image_utils::ByteMap& src,
                                    const tp_image_utils::ByteMap& labels,
                                    const CellSegmentParameters& params)
{
  size_t w = src.width();
  size_t h = src.height();

  tp_image_utils::ByteMap result = labels;

  if(result.width() != w || result.height() != h)
    return result;

  //-- Grow cells ----------------------------------------------------------------------------------
  //Here we grow the initial cells to fill all available space and try to make the most convex
  //shapes possible.
  {
    tp_image_utils::ByteMap resultB = result;
    std::array<tp_image_utils::ByteMap*, 2> buffers{&result, &resultB};
    size_t b=0;

    for(int p=0; p<params.growCellsPasses; p++)
    {
      tp_image_utils::ByteMap* srcBuffer = buffers.at(b);
      b = (b+1)%2;
      tp_image_utils::ByteMap* dstBuffer = buffers.at(b);

#if 1
      for(size_t y=0; y<h; y++)
      {
        int yoInt = (y>0)    ?-1:0;
        int yoMax = (y<(h-1))? 1:0;

        for(size_t x=0; x<w; x++)
        {
          if(src.pixel(x, y)<1)
            continue;

          int counts[8][2]={{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}};

          int xoInt = (x>0)    ?-1:0;
          int xoMax = (x<(w-1))? 1:0;
          for(int yo=yoInt; yo<=yoMax; yo++)
          {
            for(int xo=xoInt; xo<=xoMax; xo++)
            {
              if(yo==0 && xo==0)
                continue;

              int v = srcBuffer->pixel(x+size_t(xo), y+size_t(yo));

              if(v>0)
              {
                for(auto& count : counts)
                {
                  if(count[0] == 0 || count[0] == v)
                  {
                    count[0] = v;
                    count[1]++;
                    break;
                  }
                }
              }
            }
          }

          {
            int v=0;
            int c=0;
            for(auto& count : counts)
            {
              if(count[0] == 0)
                break;

              if(count[1] > c)
              {
                c = count[1];
                v = count[0];
              }
            }

            if(v>0)
              dstBuffer->setPixel(x, y, uint8_t(v));
          }
        }
      }
#endif

#if 0
      for(int y=0; y<h; y++)
      {
        int yoInt = (y>1)    ?-2:((y>0)    ?-1:0);
        int yoMax = (y<(h-2))? 2:((y<(h-1))? 1:0);

        for(int x=0; x<w; x++)
        {
          if(src.pixel(x, y)<1)
            continue;

          int counts[24][2]={{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                             {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                             {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                             {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}};

          int xoInt = (x>1)    ?-2:((x>0)    ?-1:0);
          int xoMax = (x<(w-2))? 2:((x<(w-1))? 1:0);
          for(int yo=yoInt; yo<=yoMax; yo++)
          {
            for(int xo=xoInt; xo<=xoMax; xo++)
            {
              if(yo==0 && xo==0)
                continue;

              int v = srcBuffer->pixel(x+xo, y+yo);

              if(v>0)
              {
                for(int i=0; i<24; i++)
                {
                  if(counts[i][0] == 0 || counts[i][0] == v)
                  {
                    counts[i][0] = v;
                    counts[i][1]++;
                    break;
                  }
                }
              }
            }
          }

          {
            int v=0;
            int c=0;
            for(int i=0; i<24; i++)
            {
              if(counts[i][0] == 0)
                break;

              if(counts[i][1] > c)
              {
                c = counts[i][1];
                v = counts[i][0];
              }
            }

            if(v>0)
              dstBuffer->setPixel(x, y, v);
          }
        }
      }
#endif
    }

    return *buffers.at(b);
  }
}

//##################################################################################################
tp_image_utils::ByteMap cellSegmentSimple(const tp_image_utils::ByteMap& src, const CellSegmentParameters& params)
{
  size_t w = src.width();
  size_t h = src.height();

  tp_image_utils::ByteMap result(w, h);
  result.fill(0);

  {
    //The mask will get filled out as we generate the cells
    tp_image_utils::ByteMap mask = src;

    uint8_t cellID = 0;
    uint8_t cellIDMax = uint8_t(params.maxInitialCells);
    size_t x=0;
    size_t y=0;
    size_t stride=256;

    auto innerLoops = [&]()
    {
      for(; y<h; y+=stride)
      {
        const uint8_t* r = result.constData() + (y*w);
        const uint8_t* s = src.constData()    + (y*w);

        for(x=(x<w)?x:0; x<w; x+=stride, r+=stride, s+=stride)
        {
          if((*r)==0 && (*s)==255)
          {
            cellID++;
            floodGrowCell(result, mask, cellID, int(x), int(y));
            return true;
          }
        }
      }
      return false;
    };

    for(; stride>0; stride/=2)
    {
      x=0;
      y=0;

      for(;cellID<cellIDMax;)
        if(!innerLoops())
          break;
    }
  }

  return result;
}

}
