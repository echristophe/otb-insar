/*=========================================================================

   Copyright 2009 Emmanuel Christophe
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
#include "otbForwardSensorModel.h"
#include "otbInverseSensorModel.h"
#include "otbCompositeTransform.h"
#include "otbStreamingResampleImageFilter.h"
#include "otbExtractROI.h"
#include "itkUnaryFunctorImageFilter.h"

#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

// #include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkMeanImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "otbLeeImageFilter.h"

#include "otbBinaryFunctorNeighborhoodImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "otbAmplitudePhaseToRGBFunctor.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "otbComplexToIntensityImageFilter.h"

// Command line:

// ./JustTheInterferogramComputation registered_master.tif registered_slave.tif interf.tif interf_pretty.png

/** Functor to compute the interferogram */
template< class TInput1, class TInput2, class TOutput>
    class SimpleInterferogramCalculator
{
  public:
  // The constructor and destructor.
    SimpleInterferogramCalculator() {};
    ~SimpleInterferogramCalculator() {};
  // Change detection operation
    inline TOutput operator()( const TInput1 & itA,
                            const TInput2 & itB)
    {

      TOutput result = itk::NumericTraits<TOutput>::Zero;
      double normA = 0.0;
      double normB = 0.0;
      for (unsigned long pos = 0; pos< itA.Size(); ++pos)
      {
        //TODO: add multiplication by G(m,i) function
        result += static_cast<TOutput>(itA.GetPixel(pos)* std::conj(itB.GetPixel(pos)));
        normA += itA.GetPixel(pos).real()*itA.GetPixel(pos).real()
            + itA.GetPixel(pos).imag()*itA.GetPixel(pos).imag();
        normB += itB.GetPixel(pos).real()*itB.GetPixel(pos).real()
            + itB.GetPixel(pos).imag()*itB.GetPixel(pos).imag();

      }
      if ((normA != 0) && (normB != 0))
        return static_cast<TOutput>( result/ (vcl_sqrt(normA)*vcl_sqrt(normB)) );
      else
        return itk::NumericTraits<TOutput>::Zero;
    }
};

int main(int argc, char* argv[])
{

  if (argc != 5)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile";
    std::cerr << " interferogram colorInterferogram" << std::endl;
    return EXIT_FAILURE;
    }

  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,2> ImageType;
  typedef double ScalarPixelType;
  typedef otb::Image<ScalarPixelType,2> ScalarImageType;


  /* Reading master and slave images */
  typedef otb::ImageFileReader<ImageType> ReaderType;
  typedef otb::StreamingImageFileWriter<ImageType> WriterType;

  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(argv[1]);
  slave->SetFileName(argv[2]);

  master->UpdateOutputInformation();
  slave->UpdateOutputInformation();

  std::cout << "Master size: " << std::endl;
  std::cout << master->GetOutput()->GetLargestPossibleRegion() << std::endl;

  std::cout << "Slave size: " << std::endl;
  std::cout << slave->GetOutput()->GetLargestPossibleRegion() << std::endl;

  /* Computation of the interferogram */
  typedef SimpleInterferogramCalculator<
      itk::ConstNeighborhoodIterator<ImageType>,
      itk::ConstNeighborhoodIterator<ImageType>,
      ImageType::PixelType>  InterferogramCalculatorType;

  typedef otb::BinaryFunctorNeighborhoodImageFilter<
      ImageType,ImageType,ImageType,InterferogramCalculatorType> InterferogramFilterType;
  InterferogramFilterType::Pointer interferogram = InterferogramFilterType::New();
  interferogram->SetInput1( master->GetOutput() );
  interferogram->SetInput2( slave->GetOutput() );

  ImageType::SizeType radius;
  radius[0] = 2;
  radius[1] = 3;
  interferogram->SetRadius(radius);



  WriterType::Pointer writerInterferogram = WriterType::New();
  writerInterferogram->SetFileName(argv[3]);
  writerInterferogram->SetInput(interferogram->GetOutput());
  writerInterferogram->Update();


  /* Display the interferogram with nice colors */

  
  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef otb::Image<RGBPixelType, 2> RGBImageType;

  typedef itk::ComplexToModulusImageFilter<ImageType,ScalarImageType> ModulusFilterType;
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
  writerRGB->SetFileName(argv[4]);
  writerRGB->SetInput(colormapper->GetOutput());

  writerRGB->Update();


  return EXIT_SUCCESS;
}
