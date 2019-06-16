#include "tp_image_utils_functions/SignedDistanceField.h"

#include "tp_image_utils/ColorMap.h"

#include "tp_quad_tree/QuadTreeInt.h"

#include <cmath>

namespace tp_image_utils_functions
{

namespace
{
const int cellSize = 70;
}

//##################################################################################################
tp_image_utils::ColorMap signedDistanceField(const tp_image_utils::ColorMap& src, int radius)
{
  tp_image_utils::ColorMap dst(src.width(), src.height());

  int w = src.width();
  int h = src.height();
  int maxSq = radius*radius;
  float f = 127.0f / float(radius);

  tp_quad_tree::QuadTreeInt whiteTree(0, w, 0, h, cellSize);
  tp_quad_tree::QuadTreeInt blackTree(0, w, 0, h, cellSize);

  for(int y=0; y<h; y++)
  {
    const TPPixel* s = src.constData() + (y*w);
    const TPPixel* sMax = s + (4*w);

    int x=0;
    while(s<sMax)
    {
      if(s->r>0)whiteTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
      else      blackTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
      s++;
      x++;
    }
  }

  for(int y=0; y<h; y++)
  {
    const TPPixel* s = src.constData() + (y*w);
    const TPPixel* sMax = s + w;
    TPPixel* d = dst.data() + (y*w);

    float dist=radius;

    int x=0;
    while(s<sMax)
    {
      dist++;
      int distSQ=tpMin(maxSq, int(std::ceil(dist*dist)));
      uint8_t a=0;
      if(s->r>0)
      {
        blackTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(x, y), distSQ);
        dist = std::sqrt(float(distSQ));
        a = uint8_t(127 + tpMin(dist*f, 128.0f));
      }
      else
      {
        whiteTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(x, y), distSQ);
        dist = std::sqrt(float(distSQ));
        a = uint8_t(128 - tpMin(dist*f, 128.0f));
      }

      d->r = a;
      d->g = a;
      d->b = a;
      d->a = 255;

      s++;
      d++;
      x++;
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ColorMap signedDistanceField(const tp_image_utils::ColorMap& src, int radius, int width, int height)
{
  tp_image_utils::ColorMap dst(width, height);

  int w = src.width();
  int h = src.height();
  int maxSq = radius*radius;
  float f = 127.0f / float(radius);

  tp_quad_tree::QuadTreeInt whiteTree(0, w, 0, h, cellSize);
  tp_quad_tree::QuadTreeInt blackTree(0, w, 0, h, cellSize);

  for(int y=0; y<h; y++)
  {
    const TPPixel* s = src.constData() + (y*w);
    const TPPixel* sMax = s + w;

    int x=0;
    while(s<sMax)
    {
      if(s->r>0)whiteTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
      else      blackTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
      s++;
      x++;
    }
  }

  {
    float xf = float(w)/float(width);
    float yf = float(h)/float(height);

    for(int y=0; y<height; y++)
    {
      auto sy = int(yf*float(y));
      TPPixel* d = dst.data() + (y*w);
      float dist=radius;

      for(int x=0; x<width; x++)
      {
        auto sx = int(xf*float(x));

        const TPPixel* s = src.constData() + (sy + w);
        s+=sx;

        dist+=xf;
        int distSQ=tpMin(maxSq, int(std::ceil(dist*dist)));
        uint8_t a=0;
        if(s->r>0)
        {
          blackTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(sx, sy), distSQ);
          dist = std::sqrt(float(distSQ));
          a = uint8_t(127 + tpMin(dist*f, 128.0f));
        }
        else
        {
          whiteTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(sx, sy), distSQ);
          dist = std::sqrt(float(distSQ));
          a = uint8_t(128 - tpMin(dist*f, 128.0f));
        }

        d->r = a;
        d->g = a;
        d->b = a;
        d->a = 255;

        d++;
      }
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap signedDistanceField(const tp_image_utils::ByteMap& src, int radius)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());

  int w = src.width();
  int h = src.height();
  int maxSq = radius*radius;
  float f = 127.0f / float(radius);

  tp_quad_tree::QuadTreeInt whiteTree(0, w, 0, h, cellSize);
  tp_quad_tree::QuadTreeInt blackTree(0, w, 0, h, cellSize);

  {
    const uint8_t* s = src.constData();
    for(int y=0; y<h; y++)
    {
      for(int x=0; x<w; x++)
      {
        if((*s)>0)whiteTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
        else      blackTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
        s++;
      }
    }
  }

  {
    const uint8_t* s = src.constData();
    uint8_t* d = dst.data();
    for(int y=0; y<h; y++)
    {
      float dist=radius;

      for(int x=0; x<w; x++)
      {
        dist++;
        int distSQ=tpMin(maxSq, int(std::ceil(dist*dist)));
        if((*s)>0)
        {
          blackTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(x, y), distSQ);
          dist = std::sqrt(float(distSQ));
          (*d) = uint8_t(127 + tpMin(dist*f, 128.0f));
        }
        else
        {
          whiteTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(x, y), distSQ);
          dist = std::sqrt(float(distSQ));
          (*d) = uint8_t(128 - tpMin(dist*f, 128.0f));
        }

        s++;
        d++;
      }
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap distanceField(const tp_image_utils::ByteMap& src, int radius)
{
  tp_image_utils::ByteMap dst(src.width(), src.height());

  int w = src.width();
  int h = src.height();
  int maxSq = radius*radius;
  float f = 255.0f / float(radius);

  tp_quad_tree::QuadTreeInt blackTree(0, w, 0, h, cellSize);

  {
    const uint8_t* s = src.constData();
    for(int y=0; y<h; y++)
    {
      for(int x=0; x<w; x++)
      {
        if((*s)==0)
          blackTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
        s++;
      }
    }
  }

  {
    const uint8_t* s = src.constData();
    uint8_t* d = dst.data();
    for(int y=0; y<h; y++)
    {
      //float dist=radius;

      for(int x=0; x<w; x++)
      {
        //dist++;
        //int distSQ=tpMin(maxSq, int(ceil(dist*dist)));
        if((*s)>0)
        {
          int distSQ = maxSq;
          blackTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(x, y), distSQ);
          float dist=std::sqrt(float(distSQ));
          (*d) = uint8_t(tpMin(dist*f, 255.0f));
        }
        else
        {
          (*d) = 0;
        }

        s++;
        d++;
      }
    }
  }

  return dst;
}

//##################################################################################################
tp_image_utils::ByteMap signedDistanceField(const tp_image_utils::ByteMap& src, int radius, int width, int height)
{
  tp_image_utils::ByteMap dst(width, height);

  int w = src.width();
  int h = src.height();
  int maxSq = radius*radius;
  float f = 127.0f / float(radius);

  tp_quad_tree::QuadTreeInt whiteTree(0, w, 0, h, cellSize);
  tp_quad_tree::QuadTreeInt blackTree(0, w, 0, h, cellSize);

  {
    const uint8_t* s = src.constData();
    for(int y=0; y<h; y++)
    {
      for(int x=0; x<w; x++)
      {
        if((*s)>0)whiteTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
        else      blackTree.addCoord(tp_quad_tree::QuadTreeInt::Coord(x, y));
        s++;
      }
    }
  }

  {
    float xf = float(w)/float(width);
    float yf = float(h)/float(height);

    uint8_t* d = dst.data();
    for(int y=0; y<height; y++)
    {
      auto sy = size_t(yf*float(y));
      float dist=radius;

      for(int x=0; x<width; x++)
      {
        auto sx = size_t(xf*float(x));
        dist+=xf;
        int distSQ=tpMin(maxSq, int(std::ceil(dist*dist)));
        if(src.pixel(sx, sy)>0)
        {
          blackTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(sx, sy), distSQ);
          dist = std::sqrt(float(distSQ));
          (*d) = uint8_t(127.0f + tpMin(dist*f, 128.0f));
        }
        else
        {
          whiteTree.closestPoint(tp_quad_tree::QuadTreeInt::Coord(sx, sy), distSQ);
          dist = std::sqrt(float(distSQ));
          (*d) = uint8_t(128.0f - tpMin(dist*f, 128.0f));
        }
        d++;
      }
    }
  }

  return dst;
}

}
