#ifndef tp_image_utils_functions_DeNoise_h
#define tp_image_utils_functions_DeNoise_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/ByteMap.h"

namespace tp_image_utils_functions
{

//##################################################################################################
struct ByteRegion
{
  uint8_t value{0};
  int count{0};

  int minX{0}; //!< The min x coordinate of this region calculated by calculateBoundingBoxes()
  int minY{0}; //!< The min y coordinate of this region calculated by calculateBoundingBoxes()
  int maxX{0}; //!< The max x coordinate of this region calculated by calculateBoundingBoxes()
  int maxY{0}; //!< The max y coordinate of this region calculated by calculateBoundingBoxes()
};

//##################################################################################################
struct ByteRegions
{
  //! The details for each region found
  std::vector<ByteRegion> regions;

  //! The region indexes for each pixel in the image
  std::vector<int> map;

  int w{0};
  int h{0};

  //################################################################################################
  //! Separate the regions of an image
  /*!
  Split a gray image into regions.
  \param src The image to split.
  \param addCorners Set this true if regions should be joined by corners as well as edges.
  */
  ByteRegions(const tp_image_utils::ByteMap& src, bool addCorners);

  //################################################################################################
  //! Calculate the region bounding boxes
  /*!
  Call this if you need to use use ByteRegion min and max coordinates.
  */
  void calculateBoundingBoxes();
};



//##################################################################################################
tp_image_utils::ByteMap deNoise(const tp_image_utils::ByteMap& src,
                                int minSize,
                                bool addCorners,
                                uint8_t solid=0,
                                uint8_t space=255);

//##################################################################################################
tp_image_utils::ByteMap deNoiseBlobs(const tp_image_utils::ByteMap& src,
                                     float minAspectRatio,
                                     float maxAspectRatio,
                                     float minDensity,
                                     float maxDensity,
                                     int minSize,
                                     int maxSize,
                                     bool addCorners,
                                     uint8_t solid=0,
                                     uint8_t space=255);

//##################################################################################################
tp_image_utils::ByteMap deNoiseStripes(const tp_image_utils::ByteMap& src, int minSize, uint8_t solid=0, uint8_t space=255);

//##################################################################################################
tp_image_utils::ByteMap deNoiseKnoblets(const tp_image_utils::ByteMap& src, int knobletWidth, uint8_t solid=0, uint8_t space=255);


}

#endif
