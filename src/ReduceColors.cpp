#include "tp_image_utils_functions/ReduceColors.h"

#include "tp_image_utils/ColorMap.h"

#include <unordered_map>

#include <functional>

namespace tp_image_utils_functions
{

namespace
{
//##################################################################################################
struct Color_lt
{
  int r;
  int g;
  int b;
  int count{0};

  Color_lt(int r_=0, int g_=0, int b_=0):
    r(r_),
    g(g_),
    b(b_)
  {

  }
};

bool operator==(const Color_lt& a, const Color_lt& b)
{
  return a.r==b.r && a.g==b.g && a.b==b.b;
}

//##################################################################################################
Color_lt average(const std::vector<const Color_lt*>& colors)
{
  int c=0;

  int r=0;
  int g=0;
  int b=0;

  const Color_lt*const* d = colors.data();
  const Color_lt*const* dMax = d + colors.size();

  while(d<dMax)
  {
    int w = (*d)->count;
    c+=w;
    r+=(*d)->r*w;
    g+=(*d)->g*w;
    b+=(*d)->b*w;
    d++;
  }

  if(c<1)
    c=1;

  return Color_lt(r/c, g/c, b/c);
}

//##################################################################################################
int colorDist(const Color_lt& ca, const Color_lt& cb)
{
  int r = ca.r-cb.r;
  int g = ca.g-cb.g;
  int b = ca.b-cb.b;

  return r*r + g*g + b*b;
}

//##################################################################################################
struct Bin_lt
{
  Color_lt color;
  std::vector<const Color_lt*> children;

  Bin_lt(const Color_lt& color_=Color_lt()):
    color(color_)
  {

  }
};
}
}

namespace std
{
template <>
struct hash<tp_image_utils_functions::Color_lt>
{
  std::size_t operator()(const tp_image_utils_functions::Color_lt& c) const
  {
    return std::hash<std::size_t>()(std::size_t(c.r)|(std::size_t(c.g)<<8)|(std::size_t(c.b)<<16));
  }
};
}

namespace tp_image_utils_functions
{

//##################################################################################################
tp_image_utils::ColorMap reduceColors(const tp_image_utils::ColorMap& src, int colorCount)
{
  if(colorCount<1)
    return src;

  std::unordered_map<Color_lt, int> colorHash;

  for(size_t y=0; y<src.height(); y++)
  {
    const TPPixel* s = src.constData()+(y*src.width());
    const TPPixel* sMax = s + src.width();

    while(s<sMax)
    {
      colorHash[Color_lt(s->r, s->g, s->b)]++;
      s++;
    }
  }

  //Generate an initial set of bins
  std::vector<Bin_lt> bins;
  for(int x=0; x<10; x+=4)
    for(int y=0; y<10; y+=4)
      for(int z=0; z<10; z+=4)
        bins.emplace_back(Color_lt(x*25, y*25, z*25));

  std::vector<Color_lt> colors;
  colors.reserve(colorHash.size());
  {
    for(auto i : colorHash)
    {
      Color_lt color = i.first;
      color.count = i.second;
      colors.push_back(color);
    }
  }
  const Color_lt* c = tpConst(colors).data();
  const Color_lt* cMax = c+colors.size();

  while(int(bins.size())>colorCount)
  {
    Bin_lt* binsData = bins.data();
    Bin_lt* bdMax = binsData + bins.size();

    //Empty the bins
    for(Bin_lt* bd=binsData; bd<bdMax; bd++)
      bd->children.clear();

    //Redistribute the colors
    for(const Color_lt* cc=c; cc<cMax; cc++)
    {
      int nearestDist=195075;
      Bin_lt* bdNear=binsData;
      for(Bin_lt* bd=binsData; bd<bdMax; bd++)
      {
        int dist = colorDist(bd->color, *cc);

        if(dist<nearestDist)
        {
          nearestDist=dist;
          bdNear=bd;
        }
      }

      bdNear->children.push_back(cc);
    }

    //Thin the bins
    {
      size_t least=colors.size();
      size_t removeIndex=bins.size()-1;
      size_t removeCount=0;
      for(int bb=int(removeIndex); bb>=0; bb--)
      {
        size_t count = bins.at(size_t(bb)).children.size();
        if(count==0)
        {
          bins.erase(bins.begin() + bb);
          removeCount++;
          removeIndex--;
        }
        else if(count<least)
        {
          removeIndex = size_t(bb);
          least=count;
        }
      }

      bins.erase(bins.begin() + long(removeIndex));
    }

    //Recalculate the colors
    for(int bb=0; bb<int(bins.size()); bb++)
    {
      Bin_lt& bin = bins[size_t(bb)];
      bin.color = average(bin.children);
    }
  }

  tp_image_utils::ColorMap dst(src.width(), src.height());

  {
    Bin_lt* binsData = bins.data();
    Bin_lt* bdMax = binsData + bins.size();
    const TPPixel* s = src.constData();
    const TPPixel* sMax = s + src.size();
    TPPixel* d = dst.data();

    while(s<sMax)
    {
      Color_lt original(s->r, s->g, s->b);
      Bin_lt* best=binsData;
      int nearestDist=195075;


      //Empty the bins
      for(Bin_lt* bd=binsData; bd<bdMax; bd++)
      {
        int dist = colorDist(bd->color, original);
        if(dist<nearestDist)
        {
          best=bd;
          nearestDist = dist;
        }
      }

      Color_lt& newColor = best->color;
      d->a = 255;
      d->r = uint8_t(newColor.r);
      d->g = uint8_t(newColor.g);
      d->b = uint8_t(newColor.b);

      d++;
      s++;
    }
  }

  return dst;
}

}
