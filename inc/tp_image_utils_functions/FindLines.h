#ifndef tp_image_utils_functions_FindLines_h
#define tp_image_utils_functions_FindLines_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/Point.h"
#include "tp_image_utils/ByteMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
struct FindLines
{
  //################################################################################################
  static std::vector<std::vector<tp_image_utils::Point>> findLines(const tp_image_utils::ByteMap& source,
                                                                   size_t minPoints=40,
                                                                   size_t maxDeviation=10);

  //################################################################################################
  static std::vector<std::vector<tp_image_utils::Point>> findPolylines(const tp_image_utils::ByteMap& source,
                                                                       size_t minPoints=40,
                                                                       size_t maxDeviation=10,
                                                                       size_t maxJointDistance = 100);

  //################################################################################################
  //! This returns a list of closed shapes
  /*!
  This is very similar to findPolylines, the only differences are that this will only return closed
  shapes, and the returned shape will have one fewer poins. A polyline has 1 point for each line
  plus an extra point at the end. A polygon does not need this as the last point is the same as the
  first.

  \param source - The binary image to detect the shapes in.
  \param minPoints - The minimum number of points to consider a line.
  \param maxDeviation - The max deviation from the line for a point to be considered to be part of that line.
  \param maxJointDistance - The max distance between line ends for them to be joined.

  \return A list of polygons.
  */
  static std::vector<std::vector<tp_image_utils::Point>> findPolygons(const tp_image_utils::ByteMap& source,
                                                                      size_t minPoints=40,
                                                                      size_t maxDeviation=10,
                                                                      size_t maxJointDistance = 100);

  //################################################################################################
  //! This returns a list of 4 sided closed shapes
  /*!
  This calls findPolygons and then only returns the 4 sided shapes, it also updates the point type
  to be PointTypeRectCorner.

  \param source - The binary image to detect the shapes in.
  \param minPoints - The minimum number of points to consider a line.
  \param maxDeviation - The max deviation from the line for a point to be considered to be part of that line.
  \param maxJointDistance - The max distance between line ends for them to be joined.

  \return A list of quadrilaterals.
  */
  static std::vector<std::vector<tp_image_utils::Point>> findQuadrilaterals(const tp_image_utils::ByteMap& source,
                                                                            size_t minPoints=40,
                                                                            size_t maxDeviation=10,
                                                                            size_t maxJointDistance = 100);


};

}

#endif
