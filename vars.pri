TARGET = tp_image_utils_functions
TEMPLATE = lib

DEFINES += tp_image_utils_FUNCTIONS_LIBRARY

#SOURCES += src/Globals.cpp
HEADERS += inc/tp_image_utils_functions/Globals.h


SOURCES += src/EdgeDetect.cpp
HEADERS += inc/tp_image_utils_functions/EdgeDetect.h

SOURCES += src/SignedDistanceField.cpp
HEADERS += inc/tp_image_utils_functions/SignedDistanceField.h

SOURCES += src/DeNoise.cpp
HEADERS += inc/tp_image_utils_functions/DeNoise.h

SOURCES += src/NoiseField.cpp
HEADERS += inc/tp_image_utils_functions/NoiseField.h

SOURCES += src/ConvolutionMatrix.cpp
HEADERS += inc/tp_image_utils_functions/ConvolutionMatrix.h

SOURCES += src/ReduceColors.cpp
HEADERS += inc/tp_image_utils_functions/ReduceColors.h

SOURCES += src/FindPixelGrid.cpp
HEADERS += inc/tp_image_utils_functions/FindPixelGrid.h

SOURCES += src/NormalizeBrightness.cpp
HEADERS += inc/tp_image_utils_functions/NormalizeBrightness.h

SOURCES += src/FindLines.cpp
HEADERS += inc/tp_image_utils_functions/FindLines.h

SOURCES += src/AlignImages.cpp
HEADERS += inc/tp_image_utils_functions/AlignImages.h

SOURCES += src/ExtractRect.cpp
HEADERS += inc/tp_image_utils_functions/ExtractRect.h

SOURCES += src/ExtractPolygons.cpp
HEADERS += inc/tp_image_utils_functions/ExtractPolygons.h

SOURCES += src/AddBorder.cpp
HEADERS += inc/tp_image_utils_functions/AddBorder.h

SOURCES += src/Bitwise.cpp
HEADERS += inc/tp_image_utils_functions/Bitwise.h

SOURCES += src/CellSegment.cpp
HEADERS += inc/tp_image_utils_functions/CellSegment.h

SOURCES += src/SlotFill.cpp
HEADERS += inc/tp_image_utils_functions/SlotFill.h

SOURCES += src/ToHue.cpp
HEADERS += inc/tp_image_utils_functions/ToHue.h

SOURCES += src/ToPolar.cpp
HEADERS += inc/tp_image_utils_functions/ToPolar.h

SOURCES += src/ToFloat.cpp
HEADERS += inc/tp_image_utils_functions/ToFloat.h

SOURCES += src/PixelManipulation.cpp
HEADERS += inc/tp_image_utils_functions/PixelManipulation.h

SOURCES += src/FillConcaveHull.cpp
HEADERS += inc/tp_image_utils_functions/FillConcaveHull.h

SOURCES += src/DrawLine.cpp
HEADERS += inc/tp_image_utils_functions/DrawLine.h

SOURCES += src/DrawMask.cpp
HEADERS += inc/tp_image_utils_functions/DrawMask.h

SOURCES += src/DrawShapes.cpp
HEADERS += inc/tp_image_utils_functions/DrawShapes.h

SOURCES += src/JoinImages.cpp
HEADERS += inc/tp_image_utils_functions/JoinImages.h

SOURCES += src/DilateTexture.cpp
HEADERS += inc/tp_image_utils_functions/DilateTexture.h
