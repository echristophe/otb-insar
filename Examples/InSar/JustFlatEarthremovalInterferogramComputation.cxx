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
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "otbUnaryFunctorWithIndexImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

// Command line:

// ./JustFlatEarthRemovalInterferogramComputation interf.tif phase.tif remove_earth_interf.tif cleanphase.tif


/** Functor to compute the remove earth phase */
namespace Functor
{
template< class TInput, class TOutput>
    class EarthRemovePhaseInterferogramCalculator
{
  public:
  // The constructor and destructor.
    EarthRemovePhaseInterferogramCalculator()  {};
    ~EarthRemovePhaseInterferogramCalculator() {};

	typedef itk::Index<2> IndexType;

	inline TOutput operator ()(const TInput& inPix, IndexType index)
    {
 	  TOutput result = static_cast<TOutput>( inPix);
	  for(unsigned int i = 0 ; i<2 ; ++i)
	  {
		double tmp_phase;
		tmp_phase= otb::CONST_2PI*m_IndexFrequencyMax[0]*index[0];
		TOutput phase(cos(tmp_phase), -sin(tmp_phase));
	    result *= phase;
	}
      return result;
    }

  void SetIndexFrequencyMax(IndexType value)
  {
    m_IndexFrequencyMax = value;
  }

  private:
	  IndexType m_IndexFrequencyMax;
};
}

  typedef double                           ScalarPixelType;
  typedef std::complex<ScalarPixelType>    PixelType;
  typedef otb::Image<PixelType,2>          ImageType;
  typedef ImageType::IndexType             IndexType;
  typedef otb::Image<ScalarPixelType,2>    ScalarImageType;

  typedef itk::ComplexToPhaseImageFilter<ImageType,ScalarImageType> PhaseFilterType;
  typedef itk::FFTRealToComplexConjugateImageFilter< ScalarPixelType, ImageType::ImageDimension >  RealFFTType;
  typedef itk::ComplexToModulusImageFilter<RealFFTType::OutputImageType,ScalarImageType>      ModulusFilterType;
  typedef itk::MinimumMaximumImageCalculator<ModulusFilterType::OutputImageType>  MinMaxCalculatorType;

  typedef Functor::EarthRemovePhaseInterferogramCalculator<PixelType, PixelType> InterferogramCalculatorType;
  typedef otb::UnaryFunctorWithIndexImageFilter<ImageType,ImageType,InterferogramCalculatorType> EarthRemovePhaseInterferogramFilterType;


  /* Reading images */
  typedef otb::ImageFileReader<ImageType> ReaderType;
  typedef otb::StreamingImageFileWriter<ImageType> WriterType;
  typedef otb::StreamingImageFileWriter<ScalarImageType> ScalarWriterType;

  typedef unsigned char                     OutputPixelType;
  typedef otb ::Image<OutputPixelType, 2>   OutputImageType;
  typedef itk::RescaleIntensityImageFilter<ScalarImageType , OutputImageType > RescalerType;
  typedef otb::StreamingImageFileWriter<OutputImageType > RescaleWriterType ;

int main(int argc, char* argv[])
{

  if (argc != 5)
    {
    std::cerr << "Usage: " << argv[0] << " interferogram Phase flatEarthRemoveInterferogram cleanedPhase" << std::endl;
    return EXIT_FAILURE;
    }

  ReaderType::Pointer interferogram = ReaderType::New();
  interferogram->SetFileName(argv[1]);

  PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();
  phaseFilter->SetInput(interferogram->GetOutput());
  phaseFilter->Update();

  RescalerType::Pointer rescalePhase = RescalerType::New();
  rescalePhase->SetOutputMinimum(0);
  rescalePhase->SetOutputMaximum(255);
  rescalePhase->SetInput(phaseFilter->GetOutput());

  RescaleWriterType::Pointer rescalePhaseWriter = RescaleWriterType::New();
  rescalePhaseWriter->SetFileName(argv[2]);
  rescalePhaseWriter->SetInput(rescalePhase->GetOutput());
  rescalePhaseWriter->Update();

  RealFFTType::Pointer phaseFFT = RealFFTType::New();
  phaseFFT->SetInput(phaseFilter->GetOutput());
  phaseFFT->Update();
  
  ModulusFilterType::Pointer modulusPhaseFFT = ModulusFilterType::New();
  modulusPhaseFFT->SetInput( phaseFFT->GetOutput() );
  modulusPhaseFFT->Update();

  MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
  minMax->SetImage(modulusPhaseFFT->GetOutput());
  minMax->ComputeMaximum();
  
  IndexType index;
  index = minMax->GetIndexOfMaximum();

  std::cout << "Frequency max evaluation : " << index[0]<< " with " << modulusPhaseFFT->GetOutput()->GetLargestPossibleRegion().GetSize()[0] 
	        << " , "                         << index[0]/(modulusPhaseFFT->GetOutput()->GetLargestPossibleRegion().GetSize()[0]*1.0) << std::endl;

  EarthRemovePhaseInterferogramFilterType::Pointer earthRemovePhaseInterferogram = EarthRemovePhaseInterferogramFilterType::New();
  earthRemovePhaseInterferogram->SetInput(interferogram->GetOutput());
  earthRemovePhaseInterferogram->GetFunctor().SetIndexFrequencyMax( index);
  earthRemovePhaseInterferogram->Update();

  WriterType::Pointer flatEarthRemoveWriterInterferogram = WriterType::New();
  flatEarthRemoveWriterInterferogram->SetFileName(argv[3]);
  flatEarthRemoveWriterInterferogram->SetInput(earthRemovePhaseInterferogram->GetOutput());
  flatEarthRemoveWriterInterferogram->Update();

  PhaseFilterType::Pointer cleanPhaseFilter = PhaseFilterType::New();
  cleanPhaseFilter->SetInput(earthRemovePhaseInterferogram->GetOutput());
  cleanPhaseFilter->Update();

  RescalerType::Pointer rescaleCleanPhase = RescalerType::New();
  rescaleCleanPhase->SetOutputMinimum(0);
  rescaleCleanPhase->SetOutputMaximum(255);
  rescaleCleanPhase->SetInput(cleanPhaseFilter->GetOutput());

  RescaleWriterType::Pointer rescaleCleanPhaseWriter = RescaleWriterType::New();
  rescaleCleanPhaseWriter->SetFileName(argv[4]);
  rescaleCleanPhaseWriter->SetInput(rescaleCleanPhase->GetOutput());
  rescaleCleanPhaseWriter->Update();

  return EXIT_SUCCESS;
}
