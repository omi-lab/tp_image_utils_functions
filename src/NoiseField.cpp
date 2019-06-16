#include "tp_image_utils_functions/NoiseField.h"

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ByteMap noiseField(const tp_image_utils::ByteMap& src, int radius)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());
  uint8_t* d = dst.data();
  
  int xMax = src.width();
  int yMax = src.height();
  for(int y=0; y<yMax; y++)
  {
    int startY = tpMax(y - radius, 0);
    int endY = tpMin(y + radius, yMax);
    
    for(int x=0; x<xMax; x++)
    {
      int total = 0;
      int count = 0;

      int startX = tpMax(x - radius, 0);
      int endX = tpMin(x + radius, xMax);

      for(int py=startY; py<endY; py++)
      {
        for(int px=startX; px<endX; px++)
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
      (*d) = 255-total;
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
  int xMax = src.width();
  int yMax = src.height();

  int cxMax=(xMax+(cellSize-1))/cellSize;
  int cyMax=(yMax+(cellSize-1))/cellSize;

  int binCount = cxMax*cyMax;

  std::vector<BinData_lt> bins;
  bins.resize(binCount);
  BinData_lt* binData = bins.data();

  const uint8_t* s = src.constData();

  for(int y=0; y<yMax; y++)
  {
    for(int x=0; x<xMax; x++)
    {
      BinData_lt& bin = binData[((y/cellSize)*cxMax)+(x/cellSize)];
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
    (*d) = 255-b->total;
  }

  return dst;
}
}
