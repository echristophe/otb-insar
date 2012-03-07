/*=========================================================================

   Copyright 2012 Patrick IMBO
   Contributed to ORFEO Toolbox under license Apache 2

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=========================================================================*/

// this file defines the otbCommonTest for the test driver
// and all it expects is that you have a function called RegisterTests
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "otbTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(otbBaselineFunctorBaseNew);
  REGISTER_TEST(otbHorizontalVerticalBaselineFunctorNew);
  REGISTER_TEST(otbHorizontalVerticalBaselineFunctorTest);
  REGISTER_TEST(otbLengthOrientationBaselineFunctorNew);
  REGISTER_TEST(otbLengthOrientationBaselineFunctorTest);
  REGISTER_TEST(otbParallelPerpendicularBaselineFunctorNew);
  REGISTER_TEST(otbParallelPerpendicularBaselineFunctorTest);
  REGISTER_TEST(otbPlatformPositionToBaselineCalculatorNew);
  REGISTER_TEST(otbPlatformPositionToBaselineCalculatorTest);
  REGISTER_TEST(otbBaselineCalculatorNew);
  REGISTER_TEST(otbBaselineCalculatorTest);
  REGISTER_TEST(otbMultivariateRationalTransformNew);
  REGISTER_TEST(otbMultivariateRationalTransformTest);
}
