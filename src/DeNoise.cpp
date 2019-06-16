#include "tp_image_utils_functions/DeNoise.h"

#include "tp_utils/DebugUtils.h"

#include <memory.h>

namespace tp_image_utils_functions
{

//##################################################################################################
ByteRegions::ByteRegions(const tp_image_utils::ByteMap& src, bool addCorners)
{
  w = src.width();
  h = src.height();

  if(w<1 || h<1)
    return;

  size_t ci=0;
  regions.resize(w*h);
  map.resize(w*h);
  int* done      = new int[w*h];
  auto partial = new std::pair<int, int>[w*h*8];
  size_t partialCount=0;

  memset(done, 0, w*h*sizeof(int));

  {
    const uint8_t* s = src.constData();
    const uint8_t* sMax = s + (w*h);
    int* d = map.data();
    for(; s<sMax; s++, d++)
      (*d) = (-int(*s))-1;
  }

  for(size_t y=0; y<h; y++)
  {
    for(size_t x=0; x<w; x++)
    {
      size_t offset = (y*w)+x;

      //Don't bother if we already have a count for this
      if(map[offset]>=0)
        continue;

      size_t i = ci;
      ci++;
      int color = map[offset];

      map[offset] = int(i);

      ByteRegion& region = regions[i];
      region.value=uint8_t(-1-color);
      region.count=0;

      partialCount=0;
      partial[partialCount]=std::pair<int, int>(x-1, y  ); partialCount++;
      partial[partialCount]=std::pair<int, int>(x+1, y  ); partialCount++;
      partial[partialCount]=std::pair<int, int>(x  , y-1); partialCount++;
      partial[partialCount]=std::pair<int, int>(x  , y+1); partialCount++;

      if(addCorners)
      {
        partial[partialCount]=std::pair<int, int>(x+1, y+1); partialCount++;
        partial[partialCount]=std::pair<int, int>(x+1, y-1); partialCount++;
        partial[partialCount]=std::pair<int, int>(x-1, y-1); partialCount++;
        partial[partialCount]=std::pair<int, int>(x-1, y+1); partialCount++;
      }

      while(partialCount>0)
      {
        partialCount--;
        auto px = size_t(partial[partialCount].first );
        auto py = size_t(partial[partialCount].second);

        if(px>=w || py>=h)
          continue;

        size_t po = (py*w)+px;

        if(done[po]==int(i))
          continue;

        done[po]=int(i);

        if(map[size_t(po)]!=color)
          continue;

        region.count++;
        map[po]=int(i);
        partial[partialCount]=std::pair<int, int>(px-1, py  ); partialCount++;
        partial[partialCount]=std::pair<int, int>(px+1, py  ); partialCount++;
        partial[partialCount]=std::pair<int, int>(px  , py-1); partialCount++;
        partial[partialCount]=std::pair<int, int>(px  , py+1); partialCount++;

        if(addCorners)
        {
          partial[partialCount]=std::pair<int, int>(px+1, py+1); partialCount++;
          partial[partialCount]=std::pair<int, int>(px+1, py-1); partialCount++;
          partial[partialCount]=std::pair<int, int>(px-1, py-1); partialCount++;
          partial[partialCount]=std::pair<int, int>(px-1, py+1); partialCount++;
        }
      }
    }
  }

  regions.resize(ci);

  delete[] done;
  delete[] partial;
}

//##################################################################################################
void ByteRegions::calculateBoundingBoxes()
{
  {
    ByteRegion* r = regions.data();
    ByteRegion* rMax = r + regions.size();
    for(; r<rMax; r++)
    {
      r->minX = w;
      r->minY = h;
      r->maxX = 0;
      r->maxY = 0;
    }
  }

  for(size_t y=0; y<h; y++)
  {
    for(size_t x=0; x<w; x++)
    {
      int index = map[(y*w)+x];
      ByteRegion& region = regions[size_t(index)];

      if(region.minX>x)
        region.minX=x;

      if(region.minY>y)
        region.minY=y;

      if(region.maxX<x)
        region.maxX=x;

      if(region.maxY<y)
        region.maxY=y;
    }
  }
}

//##################################################################################################
tp_image_utils::ByteMap deNoise(const tp_image_utils::ByteMap& src,
                                size_t minSize,
                                bool addCorners,
                                uint8_t solid,
                                uint8_t space)
{
  if(minSize<2)
    return src;

  size_t w = src.width();
  size_t h = src.height();

  if(w<1 || h<1)
    return src;

  tp_image_utils::ByteMap dst(w, h);
  ByteRegions regions(src, addCorners);

  uint8_t* d = dst.data();
  for(size_t y=0; y<h; y++)
  {
    for(size_t x=0; x<w; x++)
    {
      auto index = size_t(regions.map[(y*w)+x]);
      const ByteRegion& region = regions.regions[index];
      (*d) = (region.value==space || region.count<minSize)?space:solid;

      d++;
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap deNoiseBlobs(const tp_image_utils::ByteMap& src,
                     float minAspectRatio,
                     float maxAspectRatio,
                     float minDensity,
                     float maxDensity,
                     size_t minSize,
                     size_t maxSize,
                     bool addCorners,
                     uint8_t solid,
                     uint8_t space)
{
  size_t w = src.width();
  size_t h = src.height();

  if(w<1 || h<1)
    return src;

  tp_image_utils::ByteMap dst(w, h);
  ByteRegions regions(src, addCorners);
  regions.calculateBoundingBoxes();

  std::vector<int> erase;
  erase.resize(regions.regions.size());
  {
    ByteRegion* r = regions.regions.data();
    ByteRegion* rMax = r + regions.regions.size();
    int* e = erase.data();
    for(; r<rMax; r++, e++)
    {
      if(r->value==space)
      {
        (*e) = 1;
        continue;
      }

      (*e) = 0;

      size_t rw = (r->maxX - r->minX)+1;
      size_t rh = (r->maxY - r->minY)+1;

      float ar = (rw>rh)?(float(rh)/float(rw)):(float(rw)/float(rh));

      float density = float(r->count) / float(rw*rh);

      if(ar<minAspectRatio || ar>maxAspectRatio)
        continue;

      if(density<minDensity || density>maxDensity)
        continue;

      if(rw<minSize || rw>maxSize)
        continue;

      if(rh<minSize || rh>maxSize)
        continue;

      (*e) = 1;
    }
  }


  uint8_t* d = dst.data();
  for(size_t y=0; y<h; y++)
  {
    for(size_t x=0; x<w; x++)
    {
      auto index = size_t(regions.map[(y*w)+x]);
      (*d) = (erase[index])?space:solid;

      d++;
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap deNoiseStripes(const tp_image_utils::ByteMap& src,
                                       size_t minSize,
                                       uint8_t solid,
                                       uint8_t space)
{
  if(minSize<2)
    return src;

  size_t w = src.width();
  size_t h = src.height();

  if(w<1 || h<1)
    return src;

  tp_image_utils::ByteMap dst(w, h);

  //Search columns
  for(size_t x=0; x<w; x++)
  {
    size_t count=0;
    bool spaceFound=false;
    for(size_t y=0; y<h; y++)
    {
      if(src.pixel(x, y)==solid)
        count++;
      else
      {
        dst.setPixel(x, y, space);

        if(count>0)
        {
          uint8_t val = (spaceFound && count<minSize)?space:solid;
          for(size_t p=y-count; p<y; p++)
            dst.setPixel(x, p, val);

          count=0;
        }

        spaceFound=true;
      }
    }

    if(count>0)
    {
      uint8_t val = solid;
      for(size_t p=h-count; p<h; p++)
        dst.setPixel(x, p, val);
    }
  }

  //Search rows
  for(size_t y=0; y<h; y++)
  {
    size_t count=0;
    bool spaceFound=false;
    for(size_t x=0; x<w; x++)
    {
      if(dst.pixel(x, y)==solid)
        count++;
      else
      {
        dst.setPixel(x, y, space);

        if(count>0)
        {
          uint8_t val = (spaceFound && count<minSize)?space:solid;
          for(size_t p=x-count; p<x; p++)
            dst.setPixel(p, y, val);

          count=0;
        }

        spaceFound=true;
      }
    }

    if(count>0)
    {
      uint8_t val = solid;
      for(size_t p=w-count; p<w; p++)
        dst.setPixel(p, y, val);
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap deNoiseKnoblets(const tp_image_utils::ByteMap& src,
                                        size_t knobletWidth,
                                        uint8_t solid,
                                        uint8_t space)
{
  if(knobletWidth<1)
    return src;

  size_t w = src.width();
  size_t h = src.height();

  if(w<(knobletWidth+2) || h<(knobletWidth+2))
    return src;

  tp_image_utils::ByteMap dst = src;

  size_t yMax = h-(knobletWidth+1);
  for(size_t y=1; y<yMax; y++)
  {
    uint8_t* d = dst.data()+(y*w)+1;
    uint8_t* dMax = d + (w-(knobletWidth+1));
    for(; d<dMax; d++)
    {
      auto calcOnLine = [=](int xincx, int xincy, int yincy, int yincx)
      {
        auto val = [=](int x, int y)
        {
          size_t px = size_t(x*xincx) + size_t(y*yincx);
          size_t py = size_t(x*xincy) + size_t(y*yincy);

          return *(d + ((py*w) + px));
        };

        //Left side of the kernel
        if(val(-1, -1)!=space)return;
        if(val(-1,  0)!=space)return;
        if(val(-1,  1)!=solid)return;

        //Center of the kernel
        for(size_t i=0; i<=knobletWidth; i++)
        {
          if(val(int(i), -1)!=space)return;
          if(val(int(i),  1)!=solid)return;
          if(val(int(i),  0)==space)
          {
            (*d) = space;
            return;
          }
        }
      };

      auto calcOnLeftCorner = [=](int xincx, int xincy, int yincy, int yincx)
      {
        auto val = [=](int x, int y)
        {
          size_t px = size_t(x*xincx) + size_t(y*yincx);
          size_t py = size_t(x*xincy) + size_t(y*yincy);

          return *(d + (py*w) + px);
        };

        //Left side of the kernel
        if(val(-1, -1)!=space)return;
        if(val(-1,  0)!=space)return;
        if(val(-1,  1)!=space)return;

        //Center of the kernel
        for(size_t i=0; i<=knobletWidth; i++)
        {
          if(val(int(i), -1)!=space)return;
          if(val(int(i),  1)!=solid)return;
          if(val(int(i),  0)==space)
          {
            (*d) = space;
            return;
          }
        }
      };

      auto calcOnRightCorner = [=](int xincx, int xincy, int yincy, int yincx)
      {
        auto val = [=](int x, int y)
        {
          size_t px = size_t(x*xincx) + size_t(y*yincx);
          size_t py = size_t(x*xincy) + size_t(y*yincy);

          return *(d + ((py*w) + px));
        };

        //Left side of the kernel
        if(val(-1, -1)!=space)return;
        if(val(-1,  0)!=space)return;
        if(val(-1,  1)!=solid)return;

        //Center of the kernel
        for(size_t i=0; i<=knobletWidth; i++)
        {
          if(val(int(i), -1)!=space)return;
          if(val(int(i),  0)==space)
          {
            if(val(int(i),  1)==space)
              (*d) = space;
            return;
          }
          if(val(-1,  1)!=solid)return;
        }
      };

      if((*d)==space)
        continue;
      calcOnLine(1, 0,  1,  0);

      if((*d)==space)
        continue;
      calcOnLine(1, 0, -1,  0);

      if((*d)==space)
        continue;
      calcOnLine(0, 1,  0,  1);

      if((*d)==space)
        continue;
      calcOnLine(0, 1,  0, -1);

      if((*d)==space)
        continue;
      calcOnLeftCorner(1, 0,  1,  0);

      if((*d)==space)
        continue;
      calcOnLeftCorner(1, 0, -1,  0);

      if((*d)==space)
        continue;
      calcOnLeftCorner(0, 1,  0,  1);

      if((*d)==space)
        continue;
      calcOnLeftCorner(0, 1,  0, -1);

      if((*d)==space)
        continue;
      calcOnRightCorner(1, 0,  1,  0);

      if((*d)==space)
        continue;
      calcOnRightCorner(1, 0, -1,  0);

      if((*d)==space)
        continue;
      calcOnRightCorner(0, 1,  0,  1);

      if((*d)==space)
        continue;
      calcOnRightCorner(0, 1,  0, -1);
    }
  }

  return dst;
}

}
