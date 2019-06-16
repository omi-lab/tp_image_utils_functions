#include "tp_image_utils_functions/FindPixelGrid.h"

#include "tp_image_utils/Scale.h"

#include <cmath>

namespace tp_image_utils_functions
{

namespace
{
//##################################################################################################
float findBest(float calculated, float maxDiff, const std::vector<float>& points)
{
  float best = calculated;
  for(float point : points)
  {
    float diff = std::fabs(calculated - point);
    if(diff < maxDiff)
    {
      best = point;
      maxDiff = diff;
    }
  }
  return best;
}

//##################################################################################################
float calculateOffset(float cellMean, const std::vector<float>& points)
{
  std::vector<float> diffs;
  diffs.reserve(points.size());
  for(float point : points)
    diffs.emplace_back(point - (std::floor(point/cellMean)*cellMean));

  float diffsMean = 0.0f;
  float diffsSquaredSD = FindPixelGrid::squaredSD(diffs, diffsMean);
  diffsSquaredSD+=2.0f;
  return FindPixelGrid::mean(FindPixelGrid::stripOutliers(diffs, diffsMean, diffsSquaredSD));
}

//##################################################################################################
float adjustOffset(float cellMean, int cellCount, const std::vector<float>& points, float offset)
{
  if(points.empty())
    return 0.0f;

  int min = int(points.at(0) / cellMean);
  int max = int(points.at(points.size()-1) / cellMean);

  min--;
  max++;

  int best = min;
  int hits = 0;

  for(int c=min; c<=max; c++)
  {
    float minValue = (float(c)*cellMean) + offset;
    float maxValue = (float(cellCount)*cellMean) + minValue;

    minValue-=cellMean / 2.0f;
    maxValue+=cellMean / 2.0f;

    int h=0;
    for(float point : points)
      if(point>=minValue && point<=maxValue)
        h++;

    if(h>hits)
    {
      best=c;
      hits=h;
    }
  }

  return offset + (float(best)*cellMean);
}

}

//##################################################################################################
float FindPixelGrid::mean(const std::vector<int>& values)
{
  double total=0;
  for(double value : values)
    total += value;

  return float(total/double(values.size()));
}

//##################################################################################################
float FindPixelGrid::mean(const std::vector<float>& values)
{
  double total=0;
  for(float value : values)
    total += double(value);

  return float(total/double(values.size()));
}

//##################################################################################################
float FindPixelGrid::squaredSD(const std::vector<int>& values, float& mean)
{
  mean = FindPixelGrid::mean(values);

  std::vector<float> squaredDiffs;
  squaredDiffs.reserve(values.size());

  for(float value : values)
  {
    float squaredDiff = value - mean;
    squaredDiff*=squaredDiff;
    squaredDiffs.push_back(squaredDiff);
  }

  return FindPixelGrid::mean(squaredDiffs);
}

//##################################################################################################
float FindPixelGrid::squaredSD(const std::vector<float>& values, float& mean)
{
  mean = FindPixelGrid::mean(values);

  std::vector<float> squaredDiffs;
  squaredDiffs.reserve(values.size());

  for(float value : values)
  {
    float squaredDiff = value - mean;
    squaredDiff*=squaredDiff;
    squaredDiffs.push_back(squaredDiff);
  }

  return FindPixelGrid::mean(squaredDiffs);
}

//##################################################################################################
float FindPixelGrid::standardDeviation(const std::vector<int>& values, float& mean)
{
  return std::sqrt(squaredSD(values, mean));
}

//##################################################################################################
std::vector<int> FindPixelGrid::stripOutliers(const std::vector<int>& values, float mean, float squaredSD)
{
  std::vector<int> results;
  for(int value : values)
  {
    float squaredDiff = float(value) - mean;
    squaredDiff*=squaredDiff;

    if(squaredDiff<squaredSD)
      results.push_back(value);
  }

  return results;
}

//##################################################################################################
std::vector<float> FindPixelGrid::stripOutliers(const std::vector<float>& values, float mean, float squaredSD)
{
  std::vector<float> results;
  for(float value : values)
  {
    float squaredDiff = value - mean;
    squaredDiff*=squaredDiff;

    if(squaredDiff<squaredSD)
      results.push_back(value);
  }

  return results;
}

//##################################################################################################
std::vector<int> FindPixelGrid::findReversals(const std::vector<uint8_t>& src)
{
  std::vector<int> result;

  if(src.empty())
    return result;

  int8_t direction=0;
  int previous = src.at(0);
  int accumulator=0;
  bool blackSet=false;

  const uint8_t* s = src.data();
  const uint8_t* sMax = s+src.size();
  for(; s<sMax; s++)
  {
    int diff = int(*s) - previous;
    previous=(*s);

    if(diff>0)
    {
      if(direction<0)
      {
        if(blackSet)
          result.push_back(accumulator);
        accumulator=0;
        blackSet=true;
      }
    }

    accumulator++;

    if(diff)
      direction=int8_t(diff);
  }

  return result;
}

//##################################################################################################
std::vector<int> FindPixelGrid::findReversalsH(const tp_image_utils::ByteMap& src)
{
  std::vector<int> result;
  for(size_t i=0; i<src.height(); i++)
    for(int v : findReversals(src.extractColumn(i)))
      result.push_back(v);
  return result;
}

//##################################################################################################
std::vector<int> FindPixelGrid::findReversalsV(const tp_image_utils::ByteMap& src)
{
  std::vector<int> result;
  for(size_t i=0; i<src.width(); i++)
    for(int v : findReversals(src.extractColumn(i)))
      result.push_back(v);
  return result;
}

//##################################################################################################
float FindPixelGrid::findMeanPixelSize(const std::vector<int>& values)
{
  std::vector<int> reversals = values;

  float mean=0;
  float sqSD = squaredSD(reversals, mean);

  reversals = stripOutliers(reversals, mean, sqSD);

  return FindPixelGrid::mean(reversals);
}

//##################################################################################################
float FindPixelGrid::findMeanPixelWidth(const tp_image_utils::ByteMap& src)
{
  return findMeanPixelSize(findReversalsH(src));
}

//##################################################################################################
float FindPixelGrid::findMeanPixelHeight(const tp_image_utils::ByteMap& src)
{
  return findMeanPixelSize(findReversalsV(src));
}

//##################################################################################################
std::vector<uint8_t> FindPixelGrid::reversals(const std::vector<uint8_t>& src)
{
  std::vector<uint8_t> result = src;

  if(result.empty())
    return result;

  int8_t direction=0;
  int previous = result.at(0);

  uint8_t* d = result.data();
  uint8_t* dMax = d+result.size();
  for(; d<dMax; d++)
  {
    int diff = int(*d) - previous;
    previous=(*d);

    if(diff>0)
      *d = (direction<0)?0:128;
    else if(diff<0)
      *d = (direction>0)?255:128;
    else
      *d = 128;

    if(diff)
      direction=int8_t(diff);
  }

  return result;
}

//##################################################################################################
tp_image_utils::ByteMap FindPixelGrid::reversalsH(const tp_image_utils::ByteMap& src)
{
  tp_image_utils::ByteMap result = src;
  for(size_t i=0; i<src.height(); i++)
    result.setRow(i, reversals(src.extractRow(i)));
  return result;
}

//##################################################################################################
tp_image_utils::ByteMap FindPixelGrid::reversalsV(const tp_image_utils::ByteMap& src)
{
  tp_image_utils::ByteMap result = src;
  for(size_t i=0; i<src.width(); i++)
    result.setColumn(i, reversals(src.extractColumn(i)));
  return result;
}

//##################################################################################################
tp_image_utils::ByteMap FindPixelGrid::findPixelGrid(const tp_image_utils::ByteMap& src)
{
  float meanPixelWidth  = FindPixelGrid::findMeanPixelWidth(src);
  float meanPixelHeight = FindPixelGrid::findMeanPixelHeight(src);

  auto w = size_t(src.width() / meanPixelWidth);
  auto h = size_t(src.height() / meanPixelHeight);

  return scale(src, w, h);
}

//##################################################################################################
tp_image_utils::ColorMap FindPixelGrid::findPixelGrid(const tp_image_utils::ByteMap& src, const tp_image_utils::ColorMap& src2)
{
  float meanPixelWidth  = FindPixelGrid::findMeanPixelWidth(src);
  float meanPixelHeight = FindPixelGrid::findMeanPixelHeight(src);

  auto w = size_t(src.width() / meanPixelWidth);
  auto h = size_t(src.height() / meanPixelHeight);

  return scale(src2, w, h);
}

//##################################################################################################
tp_image_utils::Grid FindPixelGrid::findRegularGrid(const tp_image_utils::LineCollection& lines, tp_image_utils::GridType gridType, const FindPixelGrid::FindRegularGridParams& params)
{
  tp_image_utils::Grid grid;

  grid.type = gridType;
  grid.xCells = params.xCells;
  grid.yCells = params.yCells;

  //-- Calculate the unit vector for each line -----------------------------------------------------
  std::vector<tp_image_utils::Point> vectors;
  tp_image_utils::LineCollection actualLines;
  for(const auto& line : lines)
  {
    if(line.size()<2)
      continue;

    tp_image_utils::Point v(line.at(0).x-line.at(line.size()-1).x, line.at(0).y-line.at(line.size()-1).y);
    if(v.normalize())
    {
      vectors.push_back(v);
      actualLines.push_back(line);
    }
  }

  //-- Find the best angle -------------------------------------------------------------------------
  float bestAngle=0;
  {
    float best=0.0f;
    for(int a=0; a<90; a++)
    {
      float r = float(a)*3.1415926f/180.0f;
      tp_image_utils::Point v(std::sin(r), std::cos(r));

      float difference=0.0f;
      for(const tp_image_utils::Point& vec : vectors)
      {
        //Calculate the angle between the two vectors
        float diff = std::acos(tp_image_utils::Point::dot(vec, v));

        if(diff>1.5708f)
          diff = 3.14159f-diff;

        if(diff<0.785398f) difference += diff; // Parallel
        else difference += 1.5708f - diff;     // Perpendicular
      }

      if(difference < best || a==0)
      {
        bestAngle = a;
        best = difference;
      }
    }

    //Fine tune the angle
    float aMax = bestAngle + 1.0f;
    for(float a=bestAngle-1.0f; a<aMax; a+=0.05f)
    {
      float r = a*3.1415926f/180.0f;
      tp_image_utils::Point v(std::sin(r), std::cos(r));

      float difference=0.0f;
      for(const tp_image_utils::Point& vec : vectors)
      {
        //Calculate the angle between the two vectors
        float diff = std::acos(tp_image_utils::Point::dot(vec, v));

        if(diff>1.5708f)
          diff = 3.14159f-diff;

        if(diff<0.785398f) difference += diff; // Parallel
        else difference += 1.5708f - diff;     // Perpendicular
      }

      if(difference < best || std::fabs(a)<0.00001f)
      {
        bestAngle = a;
        best = difference;
      }
    }
  }

  //-- Filter out lines that are not parallel or perpendicular -------------------------------------
  tp_image_utils::LineCollection hLines;
  tp_image_utils::LineCollection vLines;
  {
    float t = params.angleDeviation*3.1415926f/180.0f;
    float r = float(bestAngle)*3.1415926f/180.0f;
    tp_image_utils::Point v(std::sin(r), std::cos(r));

    for(size_t i=0; i<vectors.size(); i++)
    {
      const auto& vec = vectors.at(i);

      //Calculate the angle between the two vectors
      float diff = std::acos(tp_image_utils::Point::dot(vec, v));

      if(diff>1.5708f)
        diff = 3.14159f-diff;

      if(diff<0.785398f)
      {
        if(diff<t)
          hLines.push_back(actualLines.at(i));
      }
      else
      {
        if((1.5708f - diff)<t)
          vLines.push_back(actualLines.at(i));
      }
    }
  }

  //-- Rotate the lines to orientate them with the grid --------------------------------------------
  float radians = (bestAngle-90.0f)*3.1415926f/180.0f;
  {
    rotate(hLines, radians);
    rotate(vLines, radians);
  }


  //-- Sort the lines into offsets along their axis ------------------------------------------------
  std::vector<float> hPoints; //Points along the y axis
  std::vector<float> vPoints; //Points along the x axis

  for(const tp_image_utils::Line& line : hLines)
  {
    float a=0.0f;
    for(const tp_image_utils::Point& point : line)
      a+=point.y;

    hPoints.push_back(a/float(line.size()));
  }
  for(const tp_image_utils::Line& line : vLines)
  {
    float a=0.0f;
    for(const tp_image_utils::Point& point : line)
      a+=point.x;

    vPoints.push_back(a/float(line.size()));
  }

  std::sort(hPoints.begin(), hPoints.end());
  std::sort(vPoints.begin(), vPoints.end());

  //-- Calculate the difference between each line --------------------------------------------------
  std::vector<float> hDeltas; //Deltas along the y axis
  std::vector<float> vDeltas; //Deltas along the x axis
  {
    size_t iMax = hPoints.size();
    for(size_t i=1; i<iMax; i++)
      hDeltas.push_back(hPoints.at(i)-hPoints.at(i-1));
  }
  {
    size_t iMax = vPoints.size();
    for(size_t i=1; i<iMax; i++)
      vDeltas.push_back(vPoints.at(i)-vPoints.at(i-1));
  }

  //-- Strip outliers ------------------------------------------------------------------------------
  float hMeanDelta = 0.0f;
  float vMeanDelta = 0.0f;

  float hSquaredSD = squaredSD(hDeltas, hMeanDelta);
  float vSquaredSD = squaredSD(vDeltas, vMeanDelta);

  //Add some tolerance for data sets with very little deviation. This is squared so this will have
  //very little effect of datasets with more deviation.
  hSquaredSD+=2.0f;
  vSquaredSD+=2.0f;

  std::vector<float> hDeltasFiltered = stripOutliers(hDeltas, hMeanDelta, hSquaredSD);
  std::vector<float> vDeltasFiltered = stripOutliers(vDeltas, vMeanDelta, vSquaredSD);

  if(hDeltasFiltered.empty() || vDeltasFiltered.empty())
    return grid;

  float hMean = mean(hDeltasFiltered);
  float vMean = mean(vDeltasFiltered);

  grid.xAxis.x=vMean;
  grid.xAxis.y=0.0f;

  grid.yAxis.x=0.0f;
  grid.yAxis.y=hMean;

  float hMeanHalf = hMean/2.0f;
  float vMeanHalf = vMean/2.0f;

  //-- Calculate the best offset to align this infinite grid with the lines ------------------------
  grid.origin.y = calculateOffset(hMean, hPoints);
  grid.origin.x = calculateOffset(vMean, vPoints);

  //-- Filter out lines that do not align with the grid --------------------------------------------
  {
    float xStart =  vPoints.at(0);
    float xEnd   =  vPoints.at(vPoints.size()-1);

    float yStart =  hPoints.at(0);
    float yEnd   =  hPoints.at(hPoints.size()-1);

    {
      hLines.clear();
      float threshold = hMean / 3.0f;
      for(int i=int(hPoints.size())-1; i>=0; i--)
      {
        float point = hPoints.at(size_t(i));
        float diff = point - ((std::floor(point/hMean)*hMean) + grid.origin.y);
        if(diff>hMeanHalf)
          diff-= hMean;
        if(std::fabs(diff) > threshold)
        {
          hPoints.erase(hPoints.begin()+i);
        }
        else
        {
          tp_image_utils::Line line;
          line.push_back(tp_image_utils::Point(xStart, point));
          line.push_back(tp_image_utils::Point(xEnd, point));
          hLines.push_back(line);
        }
      }
    }
    {
      vLines.clear();
      float threshold = vMean / 3.0f;
      for(int i=int(vPoints.size())-1; i>=0; i--)
      {
        float point = vPoints.at(size_t(i));
        float pOff = point - grid.origin.x;

        float diff = pOff - ((std::floor(pOff/vMean)*vMean));
        if(diff>vMeanHalf)
          diff-= vMean;

        if(std::fabs(diff) > threshold)
        {
          vPoints.erase(vPoints.begin()+i);

        }
        else
        {
          tp_image_utils::Line line;
          line.push_back(tp_image_utils::Point(point, yStart));
          line.push_back(tp_image_utils::Point(point, yEnd));
          vLines.push_back(line);
        }
      }
    }
  }

  //-- If this is a finite grid adjust that offset -------------------------------------------------
  grid.origin.y = adjustOffset(hMean, grid.yCells, hPoints, grid.origin.y);
  grid.origin.x = adjustOffset(vMean, grid.xCells, vPoints, grid.origin.x);


  //-- Calculate the distorted grid corners --------------------------------------------------------
  if(params.correctedCorners)
  {
    params.correctedCorners->clear();

    int cx = tpMax(1, grid.xCells);
    int cy = tpMax(1, grid.yCells);

    float top    = findBest(grid.origin.y,                          hMeanHalf, hPoints);
    float right  = findBest(grid.origin.x + grid.xAxis.x*float(cx), vMeanHalf, vPoints);
    float bottom = findBest(grid.origin.y + grid.yAxis.y*float(cy), hMeanHalf, hPoints);
    float left   = findBest(grid.origin.x,                          vMeanHalf, vPoints);

    params.correctedCorners->push_back(tp_image_utils::Point(left,  top   ));
    params.correctedCorners->push_back(tp_image_utils::Point(right, top   ));
    params.correctedCorners->push_back(tp_image_utils::Point(right, bottom));
    params.correctedCorners->push_back(tp_image_utils::Point(left,  bottom));
  }

  //-- Rotate the grid back to image coords --------------------------------------------------------
  grid.origin.rotate(-radians);
  grid.xAxis.rotate(-radians);
  grid.yAxis.rotate(-radians);

  if(params.hLines)
  {
    rotate(hLines, -radians);
    (*params.hLines) = hLines;
  }

  if(params.vLines)
  {
    rotate(vLines, -radians);
    (*params.vLines) = vLines;
  }

  if(params.correctedCorners)
  {
    rotate(*params.correctedCorners, -radians);

    for(auto& i : *params.correctedCorners)
      i.type = tp_image_utils::PointTypeRectCorner;
  }

  return grid;
}

}
