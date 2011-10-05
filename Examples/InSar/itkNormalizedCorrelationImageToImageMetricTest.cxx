/*=========================================================================

   Copyright 2011 Patrick IMBO
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

#include "otbImage.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

int main(int argc, char* argv[])
{


  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,2> ImageType;

  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType >		MetricType;

  MetricType::Pointer         metric        = MetricType::New();

  return EXIT_SUCCESS;
}
