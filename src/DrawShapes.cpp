#include "tp_image_utils_functions/DrawShapes.h"
#include "tp_image_utils_functions/DrawLine.h"

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ColorMap drawPoints(const tp_image_utils::ColorMap& image,
                                    const std::vector<tp_image_utils::Point>& points,
                                    const tp_image_utils::PointStyle& style)
{
  return drawPoints(image, std::vector<std::vector<tp_image_utils::Point>>({points}), style);
}

//##################################################################################################
tp_image_utils::ColorMap drawPoints(const tp_image_utils::ColorMap& image,
                                    const std::vector<std::vector<tp_image_utils::Point>>& points,
                                    const tp_image_utils::PointStyle& style)
{
  tp_image_utils::ColorMap dst = image;

  if(style.drawLine)
  {
    for(const auto& pts : points)
    {
      std::vector<std::pair<int, int>> line;

      line.reserve(pts.size()+1);
      for(const auto& p : pts)
        line.emplace_back(p.x, p.y);

      if(style.closeLine && pts.empty())
        line.emplace_back(pts.at(0).x, pts.at(0).y);

      drawPolyline(dst, style.lineColor, line);
    }
  }

  if(style.drawPoint)
  {
    for(const auto& pts : points)
    {
      for(const auto& p : pts)
      {
        int xMax = int(p.x)+3;
        int yMax = int(p.y)+3;

        for(int x=int(p.x)-2; x<xMax; x++)
          for(int y=int(p.y)-2; y<yMax; y++)
            dst.setPixel(size_t(x), size_t(y), p.color);
      }
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ColorMap drawGrid(const tp_image_utils::Grid& grid,
                                  const tp_image_utils::ColorMap& image,
                                  const tp_image_utils::GridStyle& style)
{
  int cx = grid.xCells;
  if(cx<1)
    cx=10;

  int cy = grid.yCells;
  if(cy<1)
    cy=10;

  tp_image_utils::Point xEnd(grid.xAxis.x*float(cx), grid.xAxis.y*float(cx));
  tp_image_utils::Point yEnd(grid.yAxis.x*float(cy), grid.yAxis.y*float(cy));

  //Calculate the lines that form the grid
  tp_image_utils::LineCollection lines;

  for(int c=0; c<=cx; c++)
  {
    tp_image_utils::Line line;
    line.push_back(tp_image_utils::Point(grid.origin.x + (float(c)*grid.xAxis.x),          grid.origin.y + (float(c)*grid.xAxis.y)         ));
    line.push_back(tp_image_utils::Point(grid.origin.x + yEnd.x + (float(c)*grid.xAxis.x), grid.origin.y + yEnd.y + (float(c)*grid.xAxis.y)));
    lines.push_back(line);
  }

  for(int c=0; c<=cy; c++)
  {
    tp_image_utils::Line line;
    line.push_back(tp_image_utils::Point(grid.origin.x + (float(c)*grid.yAxis.x),          grid.origin.y + (float(c)*grid.yAxis.y)         ));
    line.push_back(tp_image_utils::Point(grid.origin.x + xEnd.x + (float(c)*grid.yAxis.x), grid.origin.y + xEnd.y + (float(c)*grid.yAxis.y)));
    lines.push_back(line);
  }

  //Draw a box at the origin
  {
    tp_image_utils::Line line;
    line.push_back(tp_image_utils::Point(grid.origin.x-2, grid.origin.y-2));
    line.push_back(tp_image_utils::Point(grid.origin.x-2, grid.origin.y+2));
    line.push_back(tp_image_utils::Point(grid.origin.x+2, grid.origin.y+2));
    line.push_back(tp_image_utils::Point(grid.origin.x+2, grid.origin.y-2));
    line.push_back(tp_image_utils::Point(grid.origin.x-2, grid.origin.y-2));
    lines.push_back(line);
  }

  tp_image_utils::PointStyle ps;
  ps.lineColor = style.lineColor;
  return drawPoints(image, lines, ps);
}

}
