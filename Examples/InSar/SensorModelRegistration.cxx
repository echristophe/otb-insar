/*=========================================================================

   Copyright 2011 Julien Michel
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
#include "otbImageFileReader.h"
#include "otbStreamingImageFileWriter.h"
#include "otbGenericRSResampleImageFilter.h"
#include "otbComplexInterpolateImageFunction.h"
#include "otbWindowedSincInterpolateImageBlackmanFunction.h"
#include "otbWindowedSincInterpolateImageHammingFunction.h"
#include "otbStandardWriterWatcher.h"
#include "otbExtractROI.h"

typedef std::complex<double> PixelType;
typedef otb::Image<PixelType> ImageType;
typedef otb::ImageFileReader<ImageType> ReaderType;
typedef otb::StreamingImageFileWriter<ImageType> WriterType;
typedef otb::GenericRSResampleImageFilter<ImageType,ImageType> ResampleFilterType;
typedef otb::Function::BlackmanWindowFunction<double> FunctionType;
typedef itk::ConstantBoundaryCondition<ImageType> BoundaryConditionType;
typedef double CoordRepType;
typedef otb::ComplexInterpolateImageFunction<ImageType,FunctionType, BoundaryConditionType, CoordRepType> InterpolatorType;

typedef otb::ExtractROI<PixelType,PixelType> ExtractFilterType;


int main(int argc, char * argv[])
{
  if(argc!=10)
    {
    std::cerr<<"Usage: "<<argv[0]<<" infname reffname demdir startx starty sizex sizey outfname1 outfname2"<<std::endl;
    std::cerr<<"startx, starty, sizex and sizey defined with respect to reference image"<<std::endl;
    return EXIT_FAILURE;
    }

  ImageType::SizeType size;
  ImageType::IndexType index;

  const char * infname = argv[1];
  const char * reffname = argv[2];
  const char * demdir = argv[3];
  index[0] = atoi(argv[4]);
  index[1] = atoi(argv[5]);
  size[0] = atoi(argv[6]);
  size[1] = atoi(argv[7]);
  const char * outfname1 = argv[8];
  const char * outfname2 = argv[9];

  ReaderType::Pointer inReader = ReaderType::New();
  inReader->SetFileName(infname);
  inReader->GenerateOutputInformation();

  ReaderType::Pointer refReader = ReaderType::New();
  refReader->SetFileName(reffname);
  refReader->GenerateOutputInformation();

  ExtractFilterType::Pointer extract = ExtractFilterType::New();
  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(index);
  extract->SetExtractionRegion(region);
  extract->SetInput(refReader->GetOutput());
  extract->GetOutput()->UpdateOutputInformation();

  // Build resampler and set-up geometry
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput(inReader->GetOutput());
  resampler->SetOutputKeywordList(refReader->GetOutput()->GetImageKeywordlist());
  resampler->SetOutputOrigin(extract->GetOutput()->GetOrigin());
  resampler->SetOutputSpacing(extract->GetOutput()->GetSpacing());
  resampler->SetOutputSize(extract->GetOutput()->GetLargestPossibleRegion().GetSize());
  resampler->SetDEMDirectory(demdir);

  // Set-up interpolator
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(inReader->GetOutput());
  interpolator->SetRadius(3);
  interpolator->SetNormalizeZeroFrequency(0.01);
  resampler->SetInterpolator(interpolator);

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetInput(extract->GetOutput());
  writer1->SetFileName(outfname1);
  otb::StandardWriterWatcher watcher1(writer1,extract,"Writing extract of reference image");
  writer1->Update();
  
  
  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput(resampler->GetOutput());
  writer2->SetFileName(outfname2);

  otb::StandardWriterWatcher watcher2(writer2,resampler,"Writing registered image");

  writer2->Update();
  
  return EXIT_SUCCESS;
}
