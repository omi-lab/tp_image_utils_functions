#include "tp_image_utils_functions/DrawLine.h"

namespace tp_image_utils_functions
{

//##################################################################################################
void drawLine(tp_image_utils::ColorMap& image,
              TPPixel color,
              int x0,
              int y0,
              int x1,
              int y1)
{
  // Bresenham's line algorithm
  const bool steep = (abs(y1 - y0) > abs(x1 - x0));
  if(steep)
  {
    std::swap(x0, y0);
    std::swap(x1, y1);
  }

  if(x0 > x1)
  {
    std::swap(x0, x1);
    std::swap(y0, y1);
  }

  size_t w = image.width();
  size_t h = image.height();

  const float dx = float(x1 - x0);
  const float dy = float(std::abs(y1 - y0));

  float error = dx / 2.0f;
  const int ystep = (y0 < y1) ? 1 : -1;
  int y = y0;

  if(steep)
  {
    for(int x=x0; x<x1; x++)
    {
      if(x>=0 && y>=0 && y<int(w) && x<int(h))
        image.setPixel(size_t(y), size_t(x), color);

      error -= dy;
      if(error < 0)
      {
        y += ystep;
        error += dx;
      }
    }

  }
  else
  {
    for(int x=x0; x<x1; x++)
    {
      if(x>=0 && y>=0 && x<int(w) && y<int(h))
        image.setPixel(size_t(x), size_t(y), color);

      error -= dy;
      if(error < 0)
      {
        y += ystep;
        error += dx;
      }
    }
  }
}

//##################################################################################################
void drawPolyline(tp_image_utils::ColorMap& image,
                  TPPixel color,
                  const std::vector<std::pair<int,int>>& points)
{
  int iMax = int(points.size())-1;
  for(int i=0; i<iMax;)
  {
    const std::pair<int,int>& p0 = points.at(size_t(i));
    i++;
    const std::pair<int,int>& p1 = points.at(size_t(i));
    drawLine(image, color, p0.first, p0.second, p1.first, p1.second);
  }
}

}
