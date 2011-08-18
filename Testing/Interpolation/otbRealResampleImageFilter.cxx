/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See OTBCopyright.txt for details.


  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

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



template<class TInputPixelType, class TImageType, class TInterpolator>
int otbRealResampleImageFilter_generic(int argc, char * argv[])
{
  const char * infname = argv[1];
  const char * outfname = argv[2];

  typedef TInputPixelType											InputPixelType;
  typedef TImageType												ImageType;
  typedef TInterpolator												InterpolatorType;

  typedef InterpolatorType::ContinuousIndexType                     ContinuousIndexType;
  typedef otb::ImageFileReader<ImageType>                           ReaderType;
  typedef otb::ImageFileWriter<ImageType>                           WriterType;


  unsigned int radius = std::atoi(argv[4]);

  typedef itk::ResampleImageFilter<ImageType, ImageType, double> StreamingResampleImageFilterType;

  // Instantiating object

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(infname);
  reader->Update();

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(reader->GetOutput());
  interpolator->SetRadius(radius);
  interpolator->Initialize();

  StreamingResampleImageFilterType::Pointer resampler = StreamingResampleImageFilterType::New();
  resampler->SetInput(reader->GetOutput());
  resampler->SetInterpolator(interpolator);

  StreamingResampleImageFilterType::SizeType size;
  size[0] = 512;
  size[1] = 512;
  double tutu = 1;
  resampler->SetSize(size);
  resampler->SetOutputSpacing(tutu);

  // Result of resampler is written
  WriterType::Pointer writer     = WriterType::New();
  writer->SetInput(resampler->GetOutput());
  writer->SetFileName(outfname);
  writer->Update();


  return EXIT_SUCCESS;
}


int otbRealResampleImageFilter(int argc, char * argv[])
{
  typedef double										InputPixelType;
  typedef otb::Image<InputPixelType, 2>					ImageType;

  int FunctionType = atoi(argv[3]);
  switch (FunctionType)
    {
    case 0:
      return otbRealResampleImageFilter_generic<InputPixelType, ImageType, otb::WindowedSincInterpolateImageBlackmanFunction<ImageType> > (argc, argv);
      break;
    case 1:
      return otbRealResampleImageFilter_generic<InputPixelType, ImageType, otb::WindowedSincInterpolateImageHammingFunction<ImageType> > (argc, argv);
      break;
    default:
      std::cerr << "No more function available\n";
      return EXIT_FAILURE;
    }
}