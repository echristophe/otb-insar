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
#include "otbImageFileReader.h"
#include "otbStreamingImageFileWriter.h"
#include "otbImageFileWriter.h"
#include "itkUnaryFunctorImageFilter.h"

#include "itkTernaryFunctorImageFilter.h"
#include "otbAmplitudePhaseToRGBFunctor.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "otbComplexToIntensityImageFilter.h"
#include "itkAccumulateImageFilter.h"
#include "itkFFTComplexToComplexImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "otbFlatEarthRemovalFunctor.h"
#include "otbUnaryFunctorWithIndexImageFilter.h"

// Command line:

// ./JustFlatEarthRemovalInterferogramComputation registered_master.tif registered_slave.tif interf.tif remove_earth_interf.tif interf_pretty.png

int main(int argc, char* argv[])
{

  if (argc != 6)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile";
    std::cerr << " interferogram flatEarthRemoveInterferogram colorInterferogram" << std::endl;
    return EXIT_FAILURE;
    }

  typedef std::complex<double>    PixelType;
  typedef otb::Image<PixelType,2> ImageType;
  typedef ImageType::IndexType    IndexType;
  typedef double                  ScalarPixelType;
  typedef otb::Image<ScalarPixelType,2>                      ScalarImageType;
  typedef itk::FFTComplexToComplexImageFilter< ScalarPixelType, ImageType::ImageDimension >  FFTType;
  typedef FFTType::TransformDirectionType                                FFTDirectionType;

  typedef itk::ComplexToModulusImageFilter<FFTType::OutputImageType,ScalarImageType>      ModulusFilterType;
  typedef itk::ComplexToPhaseImageFilter<ImageType, ScalarImageType>		  FFTPhaseFilterType;
  typedef itk::MinimumMaximumImageCalculator<ScalarImageType>  MinMaxCalculatorType;

  typedef otb::Functor::FlatEarthRemovalFunctor<PixelType, PixelType> InterferogramCalculatorType;
  typedef otb::UnaryFunctorWithIndexImageFilter<ImageType,ImageType,InterferogramCalculatorType> EarthRemovePhaseInterferogramFilterType;

  /* Reading master and slave images */
  typedef otb::ImageFileReader<ImageType> ReaderType;
  typedef otb::StreamingImageFileWriter<ImageType> WriterType;

  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(argv[1]);
  slave->SetFileName(argv[2]);

  master->UpdateOutputInformation();
  slave->UpdateOutputInformation();

  ReaderType::Pointer interferogram = ReaderType::New();
  interferogram->SetFileName(argv[3]);
  interferogram->UpdateOutputInformation();

  FFTType::Pointer fft = FFTType::New();
  fft->SetInput(interferogram->GetOutput());
  fft->Update();
  
  ModulusFilterType::Pointer modulusFFT = ModulusFilterType::New();
  modulusFFT->SetInput( fft->GetOutput() );
  modulusFFT->Update();

  FFTPhaseFilterType::Pointer phaseFFT = FFTPhaseFilterType::New();
  phaseFFT->SetInput( interferogram->GetOutput() );
  phaseFFT->Update();

  MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
  minMax->SetImage(modulusFFT->GetOutput());
  minMax->ComputeMaximum();
  
  IndexType index;
  index = minMax->GetIndexOfMaximum();

  double phaseValue = phaseFFT->GetOutput()->GetPixel(index);

  std::cout << "Frequency max evaluation : " << index[0] << " , " << index[1] << std::endl;

  EarthRemovePhaseInterferogramFilterType::Pointer earthRemovePhaseInterferogram = EarthRemovePhaseInterferogramFilterType::New();
  earthRemovePhaseInterferogram->SetInput(interferogram->GetOutput());
  earthRemovePhaseInterferogram->GetFunctor().SetRangeFrequency(index[0]);
  earthRemovePhaseInterferogram->GetFunctor().SetAzimuthFrequency(index[1]);
  earthRemovePhaseInterferogram->GetFunctor().SetPhaseAtMaxFreq(phaseValue);
  earthRemovePhaseInterferogram->Update();


  WriterType::Pointer flatEarthRemoveWriterInterferogram = WriterType::New();
  flatEarthRemoveWriterInterferogram->SetFileName(argv[4]);
  flatEarthRemoveWriterInterferogram->SetInput(earthRemovePhaseInterferogram->GetOutput());
  flatEarthRemoveWriterInterferogram->Update();


  /* Display the interferogram with nice colors */

  
  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef otb::Image<RGBPixelType, 2> RGBImageType;

  ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();
  modulusFilter->SetInput(master->GetOutput());

  ModulusFilterType::Pointer coherenceFilter = ModulusFilterType::New();
  coherenceFilter->SetInput(interferogram->GetOutput());

  typedef itk::ComplexToPhaseImageFilter<ImageType,ScalarImageType> PhaseFilterType;
  PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();
  phaseFilter->SetInput(interferogram->GetOutput());

  typedef otb::Functor::AmplitudePhaseToRGBFunctor
      <ScalarPixelType,ScalarPixelType,ScalarPixelType,RGBPixelType> ColorMapFunctorType;
  typedef itk::TernaryFunctorImageFilter
      <ScalarImageType, ScalarImageType, ScalarImageType, RGBImageType, ColorMapFunctorType> ColorMapFilterType;
  ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
  colormapper->GetFunctor().SetMaximum(2000);
  colormapper->GetFunctor().SetMinimum(0);

  colormapper->SetInput1(modulusFilter->GetOutput());
  colormapper->SetInput2(coherenceFilter->GetOutput());
  colormapper->SetInput3(phaseFilter->GetOutput());
  //       colormapper->SetNumberOfThreads(1);

  typedef otb::StreamingImageFileWriter<RGBImageType> WriterRGBType;
  WriterRGBType::Pointer writerRGB = WriterRGBType::New();
  writerRGB->SetFileName(argv[5]);
  writerRGB->SetInput(colormapper->GetOutput());

  writerRGB->Update();

  return EXIT_SUCCESS;
}
