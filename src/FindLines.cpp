#include "tp_image_utils_functions/FindLines.h"

#include <unordered_set>

#include <math.h>

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
  bool taken;

  Point_lt(int x_=0, int y_=0):
    x(x_),
    y(y_),
    taken(false)
  {

  }
};

//##################################################################################################
class Histogram
{
public:
  //################################################################################################
  Histogram(int min, int max, int deviation, int maxBins):
    m_binCount(tpBound(1, (max - min)+1, maxBins)),
    m_bins(new int[m_binCount]),
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
    int* b = m_bins;
    int* bMax = b + m_binCount;
    for(; b<bMax; b++)
      (*b)=0;
  }

  //################################################################################################
  void addPoint(int value)
  {
    int i    = float((value-m_deviation) - m_min) / m_div;
    int iMax = float((value+m_deviation) - m_min) / m_div;
    iMax++;

    i    = tpBound(0, i,    m_binCount);
    iMax = tpBound(0, iMax, m_binCount);

    for(; i<iMax; i++)
      m_bins[i]++;
  }

  //################################################################################################
  //Returns the bin number
  int maxHits()
  {
    int most=0;
    int index=0;

    for(int i=0; i<m_binCount; i++)
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
  int count(int binNumber)
  {
    return (binNumber>=0 && binNumber<m_binCount)?m_bins[binNumber]:0;
  }

  //################################################################################################
  int value(int binNumber)
  {
    return (binNumber*m_div) + m_min;
  }

private:
  int   m_binCount;
  int*  m_bins;
  float m_div;
  int   m_min;
  int   m_deviation;
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
  for(const Point_lt& p : clusterPoints)
    offsetPoints.push_back(Point_lt(p.x - cx, p.y - cy));


  //-- Find the furthest point from 0 --------------------------------------------------------------
  Point_lt furthestPoint;
  float squaredDist=0;
  for(const Point_lt& p : offsetPoints)
  {
    float sq = (p.x*p.x) + (p.y*p.y);
    if(sq>squaredDist)
    {
      squaredDist = sq;
      furthestPoint = p;
    }
  }

  //-- Rotate it by 90 degrees ---------------------------------------------------------------------
  Point_lt rotatedFurthestPoint(-furthestPoint.y, furthestPoint.x);


  //-- Now rotate the points so that they are on the same side of the 90 degree line ---------------
  int iMax = offsetPoints.size();
  for(int i=0; i<iMax; i++)
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

    float length = sqrt((ax*ax) + (ay*ay));
    ax /= length;
    ay /= length;
  }

  //-- Scale the average direction -----------------------------------------------------------------
  {
    float length = sqrt((furthestPoint.x*furthestPoint.x) + (furthestPoint.y*furthestPoint.y));

    ax*=length;
    ay*=length;
  }

  std::vector<Point_lt> result;
  result.push_back(Point_lt(cx-ax, cy-ay));
  result.push_back(Point_lt(cx+ax, cy+ay));
  return result;
}

//##################################################################################################
struct IntersectionDetails_lt
{
  //An id given to each intersection so that we dont use them twice
  int id;

  //How likely is this intersection to be real
  float score;

  //The index of the lines
  int line1;
  int line2;

  //The location of the intersection
  tp_image_utils::Point point;

  IntersectionDetails_lt():
    id(-1),
    score(0.0f),
    line1(-1),
    line2(-1)
  {

  }
};

//##################################################################################################
struct LineDetails_lt
{
  std::vector<int> startIntersections;
  std::vector<int> endIntersections;

  //The index of the intersections that have been picked
  int joinStart;
  int joinEnd;

  //The extended start and end positions that we will used to calculate intersections
  tp_image_utils::Point a;
  tp_image_utils::Point b;

  bool done;

  LineDetails_lt():
    joinStart(-1),
    joinEnd(-1),
    done(false)
  {

  }
};

}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findLines(const tp_image_utils::ByteMap& source,
                                                                      int minPoints,
                                                                      int maxDeviation)
{
  std::vector<std::vector<tp_image_utils::Point> > results;

  int maxDist = (2*source.width()) + (2*source.height());

  //-- Extract the points from the source image ----------------------------------------------------
  std::vector<Point_lt> points;
  {
    const uint8_t* src = source.constData();
    int yMax = source.height();
    int xMax = source.width();
    for(int y=0; y<yMax; y++)
    {
      for(int x=0; x<xMax; x++)
      {
        if((*src)>128)
        {
          points.push_back(Point_lt(x, y));
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

  if(int(points.size())<minPoints)
    return results;

  //-- Calculate the vectors and distances ---------------------------------------------------------
  int* distances = new int[points.size() * 200];

  {
    int* dst = distances;

    //The following creates 200 vectors evenly distributed between 0 and 180 degrees. It then
    //rotates the points onto those vectors, this gives us a distance on either the x or the y axis
    //to the vector. This means that if two points have the same distance then a vector between
    //those points is parallel to the line.
    for(int c=0; c<50; c++)
    {
      float j = -(float(c+1)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = (j*float(point.y)) - float(point.x); //Distance on the x axis

      j = -(float(c)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = (j*float(point.x)) - float(point.y); //Distance on the y axis

      j = (float(c+1)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = (j*float(point.x)) - float(point.y); //Distance on the y axis

      j = (float(c)/50.0f);
      for(const Point_lt& point : tpConst(points))
        (*(dst++)) = (j*float(point.y)) - float(point.x); //Distance on the x axis
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
    int bestCount = 0;
    int bestRow   = 0;
    int bestValue = 0;

    //Find the largest cluester in the remaining data
    for(int r=0; r<200; r++)
    {
      int* dst = distances + (points.size()*r);

      int min =  maxDist;
      int max = -maxDist;
      for(int i=0; i<int(points.size()); i++)
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
      for(int i=0; i<int(points.size()); i++)
        if(!points.at(i).taken)
          histogram.addPoint(dst[i]);

      int index = histogram.maxHits();
      int count = histogram.count(index);
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
      int* dst = distances + (points.size()*bestRow);
      for(int i=0; i<int(points.size()); i++)
      {
        Point_lt& p = points[i];
        if(p.taken)
          continue;

        int deviation = abs(bestValue - dst[i]);

        if(deviation<maxDeviation)
        {
          p.taken = true;
          clusterPoints.push_back(p);
        }
      }
    }

    if(int(clusterPoints.size())<minPoints)
      break;

    //Calculate the line
    {
      std::vector<tp_image_utils::Point> line;
      for(const Point_lt& p : calculateLine(clusterPoints))
        line.push_back(tp_image_utils::Point(p.x, p.y));
      results.push_back(line);
    }
  }

  delete[] distances;
  return results;
}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findPolylines(const tp_image_utils::ByteMap& source,
                                                                          int minPoints,
                                                                          int maxDeviation,
                                                                          int maxJointDistance)
{
  float threshold = maxJointDistance*maxJointDistance;

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
  std::unordered_set<int> availableIntersections;
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
        ends.push_back({0, 1});
        ends.push_back({0, 0});
        ends.push_back({1, 0});
        ends.push_back({1, 1});

        for(int i=0; i<int(ends.size()); i++)
        {
          const std::pair<int, int>& end = ends.at(i);

          auto pLine  = line.at(end.first);
          auto pOther = other.at(end.second);

          float dx = pLine.x-pOther.x;
          float dy = pLine.y-pOther.y;

          float sqLen = (dx*dx) + (dy*dy);

          if(sqLen<threshold)
          {
            intersection.score = sqLen;
            intersection.id = intersections.size();
            availableIntersections.insert(intersection.id);
            intersections.push_back(intersection);

            ((end.first ==0)?     details.startIntersections:     details.endIntersections).push_back(intersection.id);
            ((end.second==0)?otherDetails.startIntersections:otherDetails.endIntersections).push_back(intersection.id);
          }
        }
      }
    }
  }

  //At the end of each line pick the most likely available intersection and if it is below a certain
  //threshold join those line segments. Work along the line doing this.
  while(!availableIntersections.empty())
  {
    int best=0;
    float bestScore = threshold;
    for(int index: tpConst(availableIntersections))
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
      LineDetails_lt& line1Details = lineDetails[intersection.line1];
      LineDetails_lt& line2Details = lineDetails[intersection.line2];

      if(std::find(line1Details.startIntersections.begin(),
                   line1Details.startIntersections.end(),
                   best) != line1Details.startIntersections.end())
      {
        line1Details.joinStart = best;
        for(int i : tpConst(line1Details.startIntersections))
          availableIntersections.erase(i);
      }
      else
      {
        line1Details.joinEnd = best;
        for(int i : tpConst(line1Details.endIntersections))
          availableIntersections.erase(i);
      }

      if(std::find(line2Details.startIntersections.begin(),
                   line2Details.startIntersections.end(),
                   best) != line2Details.startIntersections.end())
      {
        line2Details.joinStart = best;
        for(int i : tpConst(line2Details.startIntersections))
          availableIntersections.erase(i);
      }
      else
      {
        line2Details.joinEnd = best;
        for(int i : tpConst(line2Details.endIntersections))
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
      int idx = index;
      for(;;)
      {
        LineDetails_lt& details = lineDetails[idx];

        if(details.done)
          break;

        details.done=true;

        int join = startToEnd?details.joinEnd:details.joinStart;
        if(join>=0)
        {
          const IntersectionDetails_lt& intersection = intersections.at(join);
          result.push_back(intersection.point);
          idx = (intersection.line2==idx)?intersection.line1:intersection.line2;
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
      int idx = index;
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
          const IntersectionDetails_lt& intersection = intersections.at(join);
          result.insert(result.begin(), intersection.point);
          idx = (intersection.line2==idx)?intersection.line1:intersection.line2;
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
    if(0)
    {
      //Add the original detected lines for debug
      if(0)
        results.insert(results.end(), lines.begin(), lines.end());

      for(const IntersectionDetails_lt& intersection : tpConst(intersections))
      {
        int len=4;
        std::vector<tp_image_utils::Point> square;
        square.push_back(tp_image_utils::Point(intersection.point.x-len, intersection.point.y-len));
        square.push_back(tp_image_utils::Point(intersection.point.x-len, intersection.point.y+len));
        square.push_back(tp_image_utils::Point(intersection.point.x+len, intersection.point.y+len));
        square.push_back(tp_image_utils::Point(intersection.point.x+len, intersection.point.y-len));
        results.push_back(square);
      }
    }
  }

  return results;
}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findPolygons(const tp_image_utils::ByteMap& source,
                                                                         int minPoints,
                                                                         int maxDeviation,
                                                                         int maxJointDistance)
{
  std::vector<std::vector<tp_image_utils::Point> > polygons = FindLines::findPolylines(source, minPoints, maxDeviation, maxJointDistance);

  for(int i=polygons.size()-1; i>=0; i--)
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

    polygons.erase(polygons.begin()+i);
  }

  return polygons;
}

//##################################################################################################
std::vector<std::vector<tp_image_utils::Point> > FindLines::findQuadrilaterals(const tp_image_utils::ByteMap& source,
                                                                               int minPoints,
                                                                               int maxDeviation,
                                                                               int maxJointDistance)
{
  std::vector<std::vector<tp_image_utils::Point> > quadrilaterals = FindLines::findPolygons(source, minPoints, maxDeviation, maxJointDistance);

  for(int i=quadrilaterals.size()-1; i>=0; i--)
  {
    std::vector<tp_image_utils::Point>& quadrilateral = quadrilaterals[i];

    if(quadrilateral.size()!=4)
      quadrilaterals.erase(quadrilaterals.begin()+i);
    else
    {
      int pMax = quadrilateral.size();
      for(int p=0; p<pMax; p++)
        quadrilateral[p].type = tp_image_utils::PointTypeRectCorner;
    }
  }

  return quadrilaterals;
}

}
