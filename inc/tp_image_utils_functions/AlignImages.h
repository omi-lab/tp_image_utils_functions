#ifndef tp_image_utils_functions_AlignImages_h
#define tp_image_utils_functions_AlignImages_h

#include "tp_image_utils_functions/Globals.h"

#include "tp_image_utils/ByteMap.h"
#include "tp_image_utils/Point.h"

namespace tp_image_utils_functions
{

//##################################################################################################
struct AlignImages
{
  struct Rect
  {
    int x;
    int y;
    int w;
    int h;

    Rect(int x_=0, int y_=0, int w_=0, int h_=0):
      x(x_), y(y_), w(w_), h(h_){}

    Rect intersected(const Rect& other);
  };

  struct SkewedRegion
  {
    //! The coordinates in the reference image (the ref always remains square)
    Rect referenceRect;

    //! The 4 sided region in the other image that maps onto the referenceRect (possibly skewed)
    std::vector<tp_image_utils::Point> otherRegion;
  };

  //################################################################################################
  //! Tries to aligin two images
  /*!
  This takes two single channel images and tries to align them so the absolute subtraction of the
  two images produces the lowest result.
  */
  static std::pair<size_t, size_t> calculateMicroAlignment(const tp_image_utils::ByteMap& reference,
                                                           const tp_image_utils::ByteMap& other,
                                                           size_t maxOffset);

  //################################################################################################
  //! Translate and clip two images
  /*!
  This translates one image relative to a reference and then clips both images to only contain the
  intersection of the two images.
  */
  static void translateAndClipPair(const std::pair<int, int>& translation,
                                   tp_image_utils::ColorMap& reference,
                                   tp_image_utils::ColorMap& image);

  //################################################################################################
  //! Find the region in the other image that best fits a rect in the reference
  /*!
  This takes two single channel images and tries to align them so the absolute subtraction of the
  two images produces the lowest result.
  */
  static SkewedRegion calculateSkewedRegion(const tp_image_utils::ByteMap& reference,
                                            const tp_image_utils::ByteMap& other,
                                            size_t maxOffset);

  //################################################################################################
  //! Extract a region from one image and clip both images
  /*!
  This translates one image relative to a reference and then clips both images to only contain the
  intersection of the two images.
  */
  static void extractAndClipPair(const SkewedRegion& skewedRegion,
                                 tp_image_utils::ColorMap& reference,
                                 tp_image_utils::ColorMap& image);
};

}

#endif
