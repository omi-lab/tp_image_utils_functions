#include "tp_image_utils_functions/FindLines.h"

#include <unordered_set>

#include <cmath>

namespace tp_image_utils_functions
{

namespace
{
//##################################################################################################
struct Point_lt
{
  int x;
  int y;
  //Set true when this point has been used
  bool taken{false};

  Point_lt(int x_=0, int y_=0):
    x(x_),
    y(y_)
  {

  }
};

//##################################################################################################
class Histogram
{
  TP_NONCOPYABLE(Histogram);
public:
  //################################################################################################
  Histogram(size_t min, size_t max, size_t deviation, size_t maxBins):
    m_binCount(tpBound(size_t(1), (max - min)+1, maxBins)),
    m_bins(new size_t[m_binCount]),
    m_div(float(max-min) / float(m_binCount)),
    m_min(min),
    m_deviation(deviation)
  {
    clear();
  }

  //################################################################################################
  ~Histogram()
  {
    delete[] m_bins;
  }

  //################################################################################################
  void clear()
  {
    size_t* b = m_bins;
    size_t* bMax = b + m_binCount;
    for(; b<bMax; b++)
      (*b)=0;
  }

  //################################################################################################
  void addPoint(size_t value)
  {
    size_t i    = size_t(float((value-m_deviation) - m_min) / m_div);
    size_t iMax = size_t(float((value+m_deviation) - m_min) / m_div);
    iMax++;

    i    = tpBound(size_t(0), i,    m_binCount);
    iMax = tpBound(size_t(0), iMax, m_binCount);

    for(; i<iMax; i++)
      m_bins[i]++;
  }

  //################################################################################################
  //Returns the bin number
  size_t maxHits()
  {
    size_t most=0;
    size_t index=0;

    for(size_t i=0; i<m_binCount; i++)
    {
      if(m_bins[i]>most)
      {
        most=m_bins[i];
        index=i;
      }
    }

    return index;
  }

  //################################################################################################
  size_t count(size_t binNumber)
  {
    return (binNumber<m_binCount)?m_bins[binNumber]:0;
  }

  //################################################################################################
  size_t value(size_t binNumber)
  {
    return size_t(float(binNumber)*m_div) + m_min;
  }

private:
  size_t  m_binCount;
  size_t* m_bins;
  float   m_div;
  size_t  m_min;
  size_t  m_deviation;
};

//##################################################################################################
bool calculateIntersection(const tp_image_utils::Point& a1,
                           const tp_image_utils::Point& a2,
                           const tp_image_utils::Point& b1,
                           const tp_image_utils::Point& b2,
                           tp_image_utils::Point& i)
{
  float s1_x = a2.x - a1.x;
  float s1_y = a2.y - a1.y;

  float s2_x = b2.x - b1.x;
  float s2_y = b2.y - b1.y;

  float s = (-s1_y * (a1.x - b1.x) + s1_x * (a1.y - b1.y)) / (-s2_x * s1_y + s1_x * s2_y);
  float t = ( s2_x * (a1.y - b1.y) - s2_y * (a1.x - b1.x)) / (-s2_x * s1_y + s1_x * s2_y);

  if(s>=0 && s<=1 && t>=0 && t<=1)
  {
    i.x = a1.x + (t * s1_x);
    i.y = a1.y + (t * s1_y);
    return true;
  }

  return false;
}

//##################################################################################################
//Make a line 3 times its original length
void extendLine(tp_image_utils::Point& a, tp_image_utils::Point& b)
{
  float vx = a.x-b.x;
  float vy = a.y-b.y;
  a.x+=vx;
  a.y+=vy;
  b.x-=vx;
  b.y-=vy;
}

//##################################################################################################
std::vector<Point_lt> calculateLine(const std::vector<Point_lt>& clusterPoints)
{
  if(clusterPoints.size()<2)
    std::vector<Point_lt>();


  //-- Calculate the center point ------------------------------------------------------------------
  float cx=0;
  float cy=0;
  {
    for(const Point_lt& p : clusterPoints)
    {
      cx+=float(p.x);
      cy+=float(p.y);
    }

    cx/=float(clusterPoints.size());
    cy/=float(clusterPoints.size());
  }

  //-- Thranslate all points relative to 0 ---------------------------------------------------------
  std::vector<Point_lt> offsetPoints;
  offsetPoints.reserve(clusterPoints.size());
  for(const Point_lt& p : clusterPoints)
    offsetPoints.emplace_back(int(p.x - cx), int(p.y - cy));


  //-- Find the furthest point from 0 --------------------------------------------------------------
  Point_lt furthestPoint;
  float squaredDist=0;
  for(const Point_lt& p : offsetPoints)
  {
    float sq = float((p.x*p.x) + (p.y*p.y));
    if(sq>squaredDist)
    {
      squaredDist = sq;
      furthestPoint = p;
    }
  }

  //-- Rotate it by 90 degrees ---------------------------------------------------------------------
  Point_lt rotatedFurthestPoint(-furthestPoint.y, furthestPoint.x);


  //-- Now rotate the points so that they are on the same side of the 90 degree line ---------------
  size_t iMax = offsetPoints.size();
  for(size_t i=0; i<iMax; i++)
  {
    Point_lt& p = offsetPoints[i];
    if(((rotatedFurthestPoint.x * p.y) - (rotatedFurthestPoint.y * p.x))>0)
      p = Point_lt(-p.x, -p.y);
  }


  //-- Calculate an average direction vector -------------------------------------------------------
  float ax=0;
  float ay=0;
  {
    for(const Point_lt& p : offsetPoints)
    {
      ax+=p.x;
      ay+=p.y;
    }

    float length = std::sqrt((ax*ax) + (ay*ay));
    ax /= length;
    ay /= length;
  }

  //-- Scale the average direction -----------------------------------------------------------------
  {
    float length = std::sqrt(float((furthestPoint.x*furthestPoint.x) + (furthestPoint.y*furthestPoint.y)));

    ax*=length;
    ay*=length;
  }

  std::vector<Point_lt> result;
  result.emplace_back(int(cx-ax), int(cy-ay));
  result.emplace_back(int(cx+ax), int(cy+ay));
  return result;
}

//##################################################################################################
struct IntersectionDetails_lt
{
  //An id given to each intersection so that we dont use them twice
  int id{-1};

  //How likely is this intersection to be real
  float score{0.0f};

  //The index of the lines
  int line1{-1};
  int line2{-1};

  //The location of the intersection
  tp_image_utils::Point point;
};

//##################################################################################################
struct LineDetails_lt
{
  std::vector<size_t> startIntersections;
  std::vector<size_t> endIntersections;

  //The index of the intersections that have been picked
  int joinStart{-1};
  int joinEnd{-1};

  //The extended start and end positions that we will used to calculate intersections
  tp_image_utils::Point a;
  tp_image_utils::Point b;

  bool done{false};
};

}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findLines(const tp_image_utils::ByteMap& source,
                                                                      size_t minPoints,
                                                                      size_t maxDeviation)
{
  std::vector<std::vector<tp_image_utils::Point> > results;

  size_t maxDist = (2*source.width()) + (2*source.height());

  //-- Extract the points from the source image ----------------------------------------------------
  std::vector<Point_lt> points;
  {
    const uint8_t* src = source.constData();
    size_t yMax = source.height();
    size_t xMax = source.width();
    for(size_t y=0; y<yMax; y++)
    {
      for(size_t x=0; x<xMax; x++)
      {
        if((*src)>128)
        {
          points.emplace_back(int(x), int(y));
          if(points.size()>=10000)
          {
            y=yMax;
            break;
          }
        }

        src++;
      }
    }
  }

  if(points.size()<minPoints)
    return results;

  //-- Calculate the vectors and distances ---------------------------------------------------------
  size_t* distances = new size_t[points.size() * 200];

  {
    size_t* dst = distances;

    //The following creates 200 vectors evenly distributed between 0 and 180 degrees. It then
    //rotates the points onto those vectors, this gives us a distance on either the x or the y axis
    //to the vector. This means that if two points have the same distance then a vector between
    //those points is parallel to the line.
    for(int c=0; c<50; c++)
    {
      float j = -(float(c+1)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = size_t((j*float(point.y)) - float(point.x)); //Distance on the x axis

      j = -(float(c)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = size_t((j*float(point.x)) - float(point.y)); //Distance on the y axis

      j = (float(c+1)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = size_t((j*float(point.x)) - float(point.y)); //Distance on the y axis

      j = (float(c)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = size_t((j*float(point.y)) - float(point.x)); //Distance on the x axis
    }
  }

  //-- Find clusters -------------------------------------------------------------------------------
  //The distances array should now contain 200 rows each with points.size() values in it. What we
  //need to do now is find clusters of similar values in each row, these clusters represent points
  //that are parallel to the same vector.
  //This is done recursivly until we stop finding lines
  for(;;)
  {
    //These hold the details of the largest cluster of values
    size_t bestCount = 0;
    size_t bestRow   = 0;
    size_t bestValue = 0;

    //Find the largest cluester in the remaining data
    for(size_t r=0; r<200; r++)
    {
      size_t* dst = distances + (points.size()*r);

      size_t min = maxDist;
      size_t max = 0;
      for(size_t i=0; i<points.size(); i++)
      {
        if(points.at(i).taken)
          continue;

        if(dst[i]<min)
          min = dst[i];

        if(dst[i]>max)
          max = dst[i];
      }

      //Add all the points to a histogram
      Histogram histogram(min, max, maxDeviation/2, 10000);
      for(size_t i=0; i<points.size(); i++)
        if(!points.at(i).taken)
          histogram.addPoint(dst[i]);

      size_t index = histogram.maxHits();
      size_t count = histogram.count(index);
      if(count>bestCount)
      {
        bestCount = count;
        bestRow   = r;
        bestValue = histogram.value(index);
      }
    }

    //If count is too small then we have found all of the applicable lines
    if(bestCount<minPoints)
      break;

    //So we have a cluster of points to extract and generate a line from
    //First get the points
    std::vector<Point_lt> clusterPoints;
    {
      size_t* dst = distances + (points.size()*bestRow);
      for(size_t i=0; i<points.size(); i++)
      {
        Point_lt& p = points[i];
        if(p.taken)
          continue;

        size_t deviation = size_t(std::abs(int(bestValue) - int(dst[i])));

        if(deviation<maxDeviation)
        {
          p.taken = true;
          clusterPoints.push_back(p);
        }
      }
    }

    if(clusterPoints.size()<minPoints)
      break;

    //Calculate the line
    {
      std::vector<tp_image_utils::Point>& line = results.emplace_back();
      for(const Point_lt& p : calculateLine(clusterPoints))
        line.emplace_back(tp_image_utils::Point(float(p.x), float(p.y)));
    }
  }

  delete[] distances;
  return results;
}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findPolylines(const tp_image_utils::ByteMap& source,
                                                                          size_t minPoints,
                                                                          size_t maxDeviation,
                                                                          size_t maxJointDistance)
{
  float threshold = float(maxJointDistance*maxJointDistance);

  //-- Search for lines in the image ---------------------------------------------------------------
  std::vector<std::vector<tp_image_utils::Point>> lines = findLines(source, minPoints, maxDeviation);
  std::vector<LineDetails_lt> lineDetails;
  size_t lMax = lines.size();
  for(size_t l=0; l<lMax; l++)
  {
    const std::vector<tp_image_utils::Point>& line = lines.at(l);
    if(line.size() != 2)
    {
      lines.erase(lines.begin()+int(l));
      l--;
      lMax--;
      continue;
    }

    LineDetails_lt details;
    details.a = line.at(0);
    details.b = line.at(lines.size()-1);
    extendLine(details.a, details.b);
    lineDetails.push_back(details);
  }

  //-- Search for intersections between the lines --------------------------------------------------
  std::vector<IntersectionDetails_lt> intersections;
  std::unordered_set<size_t> availableIntersections;
  for(size_t l=0; l<lMax; l++)
  {
    auto& line = lines[l];
    LineDetails_lt& details = lineDetails[l];
    for(size_t o=l+1; o<lMax; o++)
    {
      auto& other = lines[o];
      LineDetails_lt& otherDetails = lineDetails[o];

      IntersectionDetails_lt intersection;

      if(calculateIntersection(details.a, details.b, otherDetails.a, otherDetails.b, intersection.point))
      {
        intersection.line1 = int(l);
        intersection.line2 = int(o);

        //Replicate this 4 time for (start,start)(start,end)(end,end)(end,start)
        std::vector<std::pair<int, int> > ends;
        ends.emplace_back(0, 1);
        ends.emplace_back(0, 0);
        ends.emplace_back(1, 0);
        ends.emplace_back(1, 1);

        for(const std::pair<int, int>& end : ends)
        {
          auto pLine  = line.at(size_t(end.first));
          auto pOther = other.at(size_t(end.second));

          float dx = pLine.x-pOther.x;
          float dy = pLine.y-pOther.y;

          float sqLen = (dx*dx) + (dy*dy);

          if(sqLen<threshold)
          {
            intersection.score = sqLen;
            intersection.id = int(intersections.size());
            availableIntersections.insert(size_t(intersection.id));
            intersections.push_back(intersection);

            ((end.first ==0)?     details.startIntersections:     details.endIntersections).push_back(size_t(intersection.id));
            ((end.second==0)?otherDetails.startIntersections:otherDetails.endIntersections).push_back(size_t(intersection.id));
          }
        }
      }
    }
  }

  //At the end of each line pick the most likely available intersection and if it is below a certain
  //threshold join those line segments. Work along the line doing this.
  while(!availableIntersections.empty())
  {
    size_t best=0;
    float bestScore = threshold;
    for(auto index : availableIntersections)
    {
      const IntersectionDetails_lt& intersection = intersections[index];
      if(intersection.score<bestScore)
      {
        best = index;
        bestScore = intersection.score;
      }
    }

    {
      IntersectionDetails_lt& intersection = intersections[best];
      LineDetails_lt& line1Details = lineDetails[size_t(intersection.line1)];
      LineDetails_lt& line2Details = lineDetails[size_t(intersection.line2)];

      if(std::find(line1Details.startIntersections.begin(),
                   line1Details.startIntersections.end(),
                   best) != line1Details.startIntersections.end())
      {
        line1Details.joinStart = int(best);
        for(size_t i : line1Details.startIntersections)
          availableIntersections.erase(i);
      }
      else
      {
        line1Details.joinEnd = int(best);
        for(size_t i : line1Details.endIntersections)
          availableIntersections.erase(i);
      }

      if(std::find(line2Details.startIntersections.begin(),
                   line2Details.startIntersections.end(),
                   best) != line2Details.startIntersections.end())
      {
        line2Details.joinStart = int(best);
        for(size_t i : line2Details.startIntersections)
          availableIntersections.erase(i);
      }
      else
      {
        line2Details.joinEnd = int(best);
        for(size_t i : line2Details.endIntersections)
          availableIntersections.erase(i);
      }
    }
  }

  //-- Based on the mode add polylines to the output -----------------------------------------------
  std::vector<std::vector<tp_image_utils::Point> > results;

  //Join the polylines
  for(;;)
  {
    int index = -1;
    for(size_t l=0; l<lMax; l++)
      if(!lineDetails.at(l).done)
        index = int(l);

    if(index<0)
      break;

    std::vector<tp_image_utils::Point> result;

    //Append points to this line
    {
      //Copying points from (start to end) or (end to start)
      bool startToEnd = true;
      size_t idx = size_t(index);
      for(;;)
      {
        LineDetails_lt& details = lineDetails[idx];

        if(details.done)
          break;

        details.done=true;

        int join = startToEnd?details.joinEnd:details.joinStart;
        if(join>=0)
        {
          const IntersectionDetails_lt& intersection = intersections.at(size_t(join));
          result.push_back(intersection.point);
          idx = size_t((intersection.line2==int(idx))?intersection.line1:intersection.line2);
          startToEnd = (lineDetails[idx].joinStart==join);
        }
        else
        {
          const std::vector<tp_image_utils::Point>& line = lines.at(idx);
          if(!line.empty())
            result.push_back(startToEnd?line.at(0):line.at(line.size()-1));
          break;
        }
      }
    }

    //Prepend points to this line
    {
      //Copying points from (start to end) or (end to start)
      bool startToEnd = false;
      size_t idx = size_t(index);
      lineDetails[idx].done = false;
      for(;;)
      {
        LineDetails_lt& details = lineDetails[idx];

        if(details.done)
          break;

        details.done=true;

        int join = startToEnd?details.joinEnd:details.joinStart;
        if(join>=0)
        {
          const IntersectionDetails_lt& intersection = intersections.at(size_t(join));
          result.insert(result.begin(), intersection.point);
          idx = size_t((size_t(intersection.line2)==idx)?intersection.line1:intersection.line2);
          startToEnd = (lineDetails[idx].joinEnd!=join);
        }
        else
        {
          const std::vector<tp_image_utils::Point>& line = lines.at(idx);
          result.insert(result.begin(), startToEnd?line.at(line.size()-1):line.at(0));
          break;
        }
      }
    }

    results.push_back(result);
  }

  //-- Debug results -------------------------------------------------------------------------------
  {
    //Add intersections for debug
#if 0
    {
      //Add the original detected lines for debug
#if 0
      results.insert(results.end(), lines.begin(), lines.end());
#endif
      for(const IntersectionDetails_lt& intersection : tpConst(intersections))
      {
        int len=4;
        std::vector<tp_image_utils::Point>& square = results.emplace_back();
        square.emplace_back(tp_image_utils::Point(intersection.point.x-len, intersection.point.y-len));
        square.emplace_back(tp_image_utils::Point(intersection.point.x-len, intersection.point.y+len));
        square.emplace_back(tp_image_utils::Point(intersection.point.x+len, intersection.point.y+len));
        square.emplace_back(tp_image_utils::Point(intersection.point.x+len, intersection.point.y-len));
      }
    }
#endif
  }

  return results;
}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findPolygons(const tp_image_utils::ByteMap& source,
                                                                         size_t minPoints,
                                                                         size_t maxDeviation,
                                                                         size_t maxJointDistance)
{
  std::vector<std::vector<tp_image_utils::Point> > polygons = FindLines::findPolylines(source, minPoints, maxDeviation, maxJointDistance);

  for(size_t i=polygons.size()-1; i<polygons.size(); i--)
  {
    std::vector<tp_image_utils::Point>& polygon = polygons[i];

    if(polygon.size()>3)
    {
      const tp_image_utils::Point& start = polygon.at(0);
      const tp_image_utils::Point& end   = polygon.at(polygon.size()-1);

      float dx = start.x - end.x;
      float dy = start.y - end.y;
      float dist = (dx*dx) + (dy*dy);
      if(dist <1.5f)
      {
        polygon.erase(polygon.end()-1);
        continue;
      }
    }

    polygons.erase(polygons.begin()+ptrdiff_t(i));
  }

  return polygons;
}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findQuadrilaterals(const tp_image_utils::ByteMap& source,
                                                                               size_t minPoints,
                                                                               size_t maxDeviation,
                                                                               size_t maxJointDistance)
{
  std::vector<std::vector<tp_image_utils::Point> > quadrilaterals = FindLines::findPolygons(source, minPoints, maxDeviation, maxJointDistance);

  for(size_t i=quadrilaterals.size()-1; i<quadrilaterals.size(); i--)
  {
    std::vector<tp_image_utils::Point>& quadrilateral = quadrilaterals[i];

    if(quadrilateral.size()!=4)
      quadrilaterals.erase(quadrilaterals.begin()+ptrdiff_t(i));
    else
    {
      size_t pMax = quadrilateral.size();
      for(size_t p=0; p<pMax; p++)
        quadrilateral[p].type = tp_image_utils::PointTypeRectCorner;
    }
  }

  return quadrilaterals;
}

}
