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
#include "itkMacro.h"

#include "otbComplexInterpolateImageFunction.h"
#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "otbStreamingImageFileWriter.h"
#include <fstream>

#include "otbWindowedSincInterpolateImageBlackmanFunction.h"
#include "otbWindowedSincInterpolateImageHammingFunction.h"

template<class TFunction>
int otbComplexResampleImageFilter_generic(int argc, char * argv[])
{
  const char * infname = argv[1];
  const char * outfname = argv[2];

  typedef std::complex<double>                    InputPixelType;
  typedef otb::Image<InputPixelType, 2>                ImageType;
  typedef TFunction                            FunctionType;
  typedef itk::ConstantBoundaryCondition<ImageType>          BoundaryConditionType;
  typedef double                          CoordRepType;

  typedef otb::ComplexInterpolateImageFunction<ImageType, 
            FunctionType, BoundaryConditionType,
            CoordRepType>                InterpolatorType;

  typedef typename InterpolatorType::ContinuousIndexType            ContinuousIndexType;
  typedef otb::ImageFileReader<ImageType>                           ReaderType;
  typedef otb::ImageFileWriter<ImageType>                           WriterType;


  unsigned int radius = std::atoi(argv[4]);
  double zeroFrequency = std::atof(argv[5]);

  typedef itk::ResampleImageFilter<ImageType, ImageType, double> StreamingResampleImageFilterType;

  // Instantiating object

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(infname);
  reader->Update();

  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(reader->GetOutput());
  interpolator->SetRadius(radius);
  interpolator->SetNormalizeZeroFrequency(zeroFrequency);

  typename StreamingResampleImageFilterType::Pointer resampler = StreamingResampleImageFilterType::New();
  resampler->SetInput(reader->GetOutput());
  resampler->SetInterpolator(interpolator);

  typename StreamingResampleImageFilterType::SizeType size;
  size[0] = 512;
  size[1] = 512;
  double tutu = 1;
  resampler->SetSize(size);
  resampler->SetOutputSpacing(tutu);

  // Result of resampler is written
  typename WriterType::Pointer writer     = WriterType::New();
  writer->SetInput(resampler->GetOutput());
  writer->SetFileName(outfname);
  writer->Update();

  return EXIT_SUCCESS;
}


int otbComplexResampleImageFilter(int argc, char * argv[])
{
  typedef std::complex<double>                    InputPixelType;
  typedef otb::Image<InputPixelType, 2>                ImageType;

  int FunctionType = atoi(argv[3]);
  switch (FunctionType)
    {
    case 0:
      return otbComplexResampleImageFilter_generic<otb::Function::BlackmanWindowFunction<double> > (argc, argv);
      break;
    case 1:
      return otbComplexResampleImageFilter_generic<otb::Function::HammingWindowFunction<double> > (argc, argv);
      break;
    default:
      std::cerr << "No more function available\n";
      return EXIT_FAILURE;
    }
}
