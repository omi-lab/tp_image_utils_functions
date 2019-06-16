#include "tp_image_utils_functions/ToFloat.h"

namespace tp_image_utils_functions
{

//##################################################################################################
std::vector<std::string> channelModes()
{
  return {"Interleaved", "Separate"};
}

//##################################################################################################
std::string channelModeToString(ChannelMode channelMode)
{
  switch(channelMode)
  {
  case ChannelMode::Interleaved: return "Interleaved";
  case ChannelMode::Separate:    return "Separate";
  }
  return "Interleaved";
}

//##################################################################################################
ChannelMode channelModeFromString(const std::string& channelMode)
{
  if(channelMode == "Interleaved") return ChannelMode::Interleaved;
  return ChannelMode::Separate;
}

//##################################################################################################
std::vector<std::string> channelOrders()
{
  return {"RGB", "BGR"};
}

//##################################################################################################
std::string channelOrderToString(ChannelOrder channelOrder)
{
  switch(channelOrder)
  {
  case ChannelOrder::RGB: return "RGB";
  case ChannelOrder::BGR: return "BGR";
  }
  return "RGB";
}

//##################################################################################################
ChannelOrder channelOrderFromString(const std::string& channelOrder)
{
  if(channelOrder == "RGB") return ChannelOrder::RGB;
  return ChannelOrder::BGR;
}

//##################################################################################################
void toFloat(const tp_image_utils::ColorMap& src,
           ChannelMode channelMode,
           ChannelOrder channelOrder,
           std::vector<float>& outData)
{
  size_t size = src.width() * src.height();
  outData.resize(size*3);

  auto i = src.constData();
  auto iMax = i+size;

  switch (channelMode)
  {
  case tp_image_utils_functions::ChannelMode::Interleaved:
  {
    switch(channelOrder)
    {
    case tp_image_utils_functions::ChannelOrder::RGB:
    {
      for(auto j = outData.data(); i<iMax; i++, j+=3)
      {
        j[0] = float(i->r) / 255.0f;
        j[1] = float(i->g) / 255.0f;
        j[2] = float(i->b) / 255.0f;
      }
      break;
    }

    case tp_image_utils_functions::ChannelOrder::BGR:
    {
      for(auto j = outData.data(); i<iMax; i++, j+=3)
      {
        j[0] = float(i->b) / 255.0f;
        j[1] = float(i->g) / 255.0f;
        j[2] = float(i->r) / 255.0f;
      }
      break;
    }
    }
    break;
  }

  case tp_image_utils_functions::ChannelMode::Separate:
  {
    float* ro{nullptr};
    float* go{nullptr};
    float* bo{nullptr};

    switch(channelOrder)
    {
    case tp_image_utils_functions::ChannelOrder::RGB:
    {
      ro = outData.data();
      go = ro + size;
      bo = go + size;
      break;
    }

    case tp_image_utils_functions::ChannelOrder::BGR:
    {
      bo = outData.data();
      go = bo + size;
      ro = go + size;
      break;
    }
    }

    for(; i<iMax; i++, ro++, go++, bo++)
    {
      (*ro) = float(i->r) / 255.0f;
      (*go) = float(i->g) / 255.0f;
      (*bo) = float(i->b) / 255.0f;
    }
    break;
  }
  }
}

}
