#include "tp_image_utils_functions/PixelManipulation.h"

#if defined(TP_WIN32_MINGW) && defined(TP_DEBUG)
#define NO_EXPRTK
#endif

#ifndef NO_EXPRTK
#include "exprtk.hpp"
#endif

namespace tp_image_utils_functions
{

namespace
{
#ifndef NO_EXPRTK
//##################################################################################################
using symbol_table_lt =  exprtk::symbol_table<float>;
using   expression_lt =  exprtk::expression<float>  ;
using       parser_lt =  exprtk::parser<float>      ;

//##################################################################################################
struct Expr_lt
{
  symbol_table_lt symbolTable;
  expression_lt expression;

  expression_lt calcRed;
  expression_lt calcGreen;
  expression_lt calcBlue;
  expression_lt calcAlpha;
  expression_lt calcByte;

  float red   {1.0f};
  float green {1.0f};
  float blue  {1.0f};
  float alpha {1.0f};
  float byte  {1.0f};
};

//##################################################################################################
void populateIn(const TPPixel* in, Expr_lt& e)
{
  e.red   = float(in->r) / 255.0f;
  e.green = float(in->g) / 255.0f;
  e.blue  = float(in->b) / 255.0f;
  e.alpha = float(in->a) / 255.0f;
}

//##################################################################################################
void populateIn(const uint8_t* in, Expr_lt& e)
{
  e.byte = float(*in) / 255.0f;
}

//##################################################################################################
void calculateOut(TPPixel* out, Expr_lt& e)
{
  out->r = uint8_t(e.calcRed  .value()*255.0f);
  out->g = uint8_t(e.calcGreen.value()*255.0f);
  out->b = uint8_t(e.calcBlue .value()*255.0f);
  out->a = uint8_t(e.calcAlpha.value()*255.0f);
}

//##################################################################################################
void calculateOut(uint8_t* out, Expr_lt& e)
{
  (*out) = uint8_t(e.calcByte.value()*255.0f);
}

//##################################################################################################
void parse(const std::string& expressionString, expression_lt& expression, std::vector<std::string>& errors)
{
  parser_lt parser;
  if(!parser.compile(expressionString,expression))
  {
    errors.push_back("Error: " + parser.error() + " Expression: " + expressionString);

    for (size_t i = 0; i < parser.error_count(); ++i)
    {
      const auto error = parser.get_error(i);
      errors.push_back("Error: " + std::to_string(static_cast<int>(i)) +
                       " Position: " + std::to_string(static_cast<int>(error.token.position)) +
                       " Type: " + exprtk::parser_error::to_str(error.mode) +
                       " Msg: " + error.diagnostic +
                       " Expr: " + expressionString);
    }
  }
}

//##################################################################################################
template<typename Out, typename In>
Out pixelManipulation(const In& in, const PixelManipulation& params, std::vector<std::string>& errors)
{
  Out out(in.width(), in.height());

  Expr_lt e;

  e.symbolTable.add_constants();
  e.symbolTable.add_variable("red"  , e.red  );
  e.symbolTable.add_variable("green", e.green);
  e.symbolTable.add_variable("blue" , e.blue );
  e.symbolTable.add_variable("alpha", e.alpha);
  e.symbolTable.add_variable("byte" , e.byte );

  e.calcRed  .register_symbol_table(e.symbolTable);
  e.calcGreen.register_symbol_table(e.symbolTable);
  e.calcBlue .register_symbol_table(e.symbolTable);
  e.calcAlpha.register_symbol_table(e.symbolTable);
  e.calcByte .register_symbol_table(e.symbolTable);

  parse(params.calcRed  , e.calcRed  , errors);
  parse(params.calcGreen, e.calcGreen, errors);
  parse(params.calcBlue , e.calcBlue , errors);
  parse(params.calcAlpha, e.calcAlpha, errors);
  parse(params.calcByte , e.calcByte , errors);

  if(errors.empty())
  {
    auto src = in.constData();
    auto srcMax = src+in.size();
    auto dst = out.data();
    for(; src<srcMax; src++, dst++)
    {
      populateIn(src, e);
      calculateOut(dst, e);
    }
  }

  return out;
}
#else
//##################################################################################################
// This won't compile on MinGW debug builds so just stubb it.
template<typename Out, typename In>
Out pixelManipulation(const In& in, const PixelManipulation& params, std::vector<std::string>& errors)
{
  TP_UNUSED(params);
  TP_UNUSED(errors);
  Out out(in.width(), in.height());
  return out;
}
#endif
}

//##################################################################################################
tp_image_utils::ColorMap pixelManipulationColor(const tp_image_utils::ColorMap& src, const PixelManipulation& params, std::vector<std::string>& errors)
{
  return pixelManipulation<tp_image_utils::ColorMap>(src, params, errors);
}

//##################################################################################################
tp_image_utils::ColorMap pixelManipulationColor(const tp_image_utils::ByteMap& src, const PixelManipulation& params, std::vector<std::string>& errors)
{
  return pixelManipulation<tp_image_utils::ColorMap>(src, params, errors);
}

//##################################################################################################
tp_image_utils::ByteMap pixelManipulationByte(const tp_image_utils::ColorMap& src, const PixelManipulation& params, std::vector<std::string>& errors)
{
  return pixelManipulation<tp_image_utils::ByteMap>(src, params, errors);
}

//##################################################################################################
tp_image_utils::ByteMap pixelManipulationByte(const tp_image_utils::ByteMap& src, const PixelManipulation& params, std::vector<std::string>& errors)
{
  return pixelManipulation<tp_image_utils::ByteMap>(src, params, errors);
}

}
