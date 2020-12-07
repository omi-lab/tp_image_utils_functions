#include "tp_image_utils_functions/NoiseField.h"

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ByteMap noiseField(const tp_image_utils::ByteMap& src, int radius)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());
  uint8_t* d = dst.data();
  
  size_t xMax = src.width();
  size_t yMax = src.height();
  for(size_t y=0; y<yMax; y++)
  {
    size_t startY = tpMax(y - size_t(radius), size_t(0));
    size_t endY = tpMin(y + size_t(radius), yMax);
    
    for(size_t x=0; x<xMax; x++)
    {
      int total = 0;
      int count = 0;

      size_t startX = tpMax(x - size_t(radius), size_t(0));
      size_t endX = tpMin(x + size_t(radius), xMax);

      for(size_t py=startY; py<endY; py++)
      {
        for(size_t px=startX; px<endX; px++)
        {
          if(src.pixel(px, py)>0)
            total++;
          count++;
        }
      }

      int cH = count/2;
      if(cH<1)
        cH=1;

      total = abs(total-cH);
      total = (total*255)/cH;
      (*d) = uint8_t(255-total);
      d++;
    }
  }

  return dst;
}

namespace
{
struct BinData_lt
{
  int count{0};
  int total{0};
};
}

//##################################################################################################
tp_image_utils::ByteMap noiseFieldGrid(const tp_image_utils::ByteMap& src, int cellSize)
{
  size_t xMax = src.width();
  size_t yMax = src.height();

  size_t cxMax=(xMax+tpMax(size_t(1), size_t(cellSize))-size_t(1))/size_t(cellSize);
  size_t cyMax=(yMax+tpMax(size_t(1), size_t(cellSize))-size_t(1))/size_t(cellSize);

  size_t binCount = cxMax*cyMax;

  std::vector<BinData_lt> bins;
  bins.resize(binCount);
  BinData_lt* binData = bins.data();

  const uint8_t* s = src.constData();

  for(size_t y=0; y<yMax; y++)
  {
    for(size_t x=0; x<xMax; x++)
    {
      BinData_lt& bin = binData[((y/size_t(cellSize))*cxMax)+(x/size_t(cellSize))];
      if((*s)>0)
        bin.total++;
      bin.count++;
      s++;
    }
  }

  tp_image_utils::ByteMap dst(cxMax, cyMax);
  uint8_t* d = dst.data();
  BinData_lt* b = binData;
  uint8_t* dMax = d + binCount;
  for(; d<dMax; d++, b++)
  {
    int cH = b->count/2;

    b->total = abs(b->total-cH);
    b->total = (b->total*255)/cH;
    (*d) = uint8_t(255-b->total);
  }

  return dst;
}
}
