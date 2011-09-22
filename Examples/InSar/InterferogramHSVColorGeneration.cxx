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
#include "itkTernaryFunctorImageFilter.h"

#include "itkRGBPixel.h"
#include "otbScalarToRainbowRGBPixelFunctor.h"
#include "otbMath.h"

#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkMinimumMaximumImageCalculator.h"


// Command line:

// ./JustTheInterferogramComputation registered_master.tif registered_slave.tif interf.tif interf_pretty.png


template<class TInput1, class TInput2 = TInput1, class TInput3 = TInput1, class TOutput = TInput1>
class AmplitudeCoherencyPhaseToRGBFunctor
{
public:
  typedef TOutput                          RGBPixelType;
  typedef typename RGBPixelType::ValueType RGBComponentType;
  typedef otb::Functor::HSVToRGBFunctor<RGBPixelType>    HSVToRGBFunctorType;
  typedef TInput1                          ScalarType;

  AmplitudeCoherencyPhaseToRGBFunctor()
    {
    m_Minimum = 0;
    m_Maximum = itk::NumericTraits<ScalarType>::max();
    };
  ~AmplitudeCoherencyPhaseToRGBFunctor(){}

  void SetMaximum(ScalarType max)
  {
    this->m_Maximum = max;
  }

  void SetMinimum(ScalarType min)
  {
    this->m_Minimum = min;
  }

  inline TOutput operator ()(const TInput1& amplitude, const TInput2& coherence, const TInput3& phase) const
  {
    double hinc, sinc, vinc;
    hinc = 1.0	 / (otb::CONST_2PI);
    sinc = 0.0;
    vinc = 0.0;

    double hue, sat, val;

    hue = (phase + otb::CONST_PI) * hinc;
    sat = coherence;
    val = itk::NumericTraits<RGBComponentType>::max() / 2.
          * ((amplitude - m_Minimum) / (m_Maximum - m_Minimum));

    if (amplitude < m_Minimum)
      {
      val = 0;
      }
    if (amplitude > m_Maximum)
      {
      val = itk::NumericTraits<RGBComponentType>::max();
      }

    return m_HSVToRGBFunctor(hue, sat, val);

  }
private:
  ScalarType          m_Maximum;
  ScalarType          m_Minimum;
  HSVToRGBFunctorType m_HSVToRGBFunctor;
};




int main(int argc, char* argv[])
{

  if (argc != 5)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile";
    std::cerr << " interferogram HSVcolorImage" << std::endl;
    return EXIT_FAILURE;
    }

  typedef double                           ScalarPixelType;
  typedef std::complex<ScalarPixelType>    PixelType;
  typedef otb::Image<PixelType,2>          ImageType;
  typedef ImageType::IndexType             IndexType;
  typedef ImageType::SizeType              SizeType;
  typedef otb::Image<ScalarPixelType,2>    ScalarImageType;

  /* Reading master and slave images */
  typedef otb::ImageFileReader<ImageType> ReaderType;
  typedef otb::StreamingImageFileWriter<ImageType> WriterType;

  typedef itk::ComplexToModulusImageFilter<ImageType,ScalarImageType> ModulusFilterType;
  typedef itk::Statistics::ScalarImageToHistogramGenerator<ScalarImageType> HistogramGeneratorType;

  ReaderType::Pointer master = ReaderType::New();
  master->SetFileName(argv[1]);
  master->UpdateOutputInformation();

  ReaderType::Pointer slave = ReaderType::New();
  slave->SetFileName(argv[2]);
  slave->UpdateOutputInformation();

  ReaderType::Pointer interferogram = ReaderType::New();
  interferogram->SetFileName(argv[3]);
  interferogram->UpdateOutputInformation();

  ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();
  modulusFilter->SetInput(master->GetOutput());
  modulusFilter->Update();

  HistogramGeneratorType::Pointer modulusHistogram = HistogramGeneratorType::New();
  modulusHistogram->SetInput(modulusFilter->GetOutput());
  modulusHistogram->SetNumberOfBins( 500.0 );
  modulusHistogram->SetMarginalScale( 10.0 );
  modulusHistogram->Compute();

  double minModulus = modulusHistogram->GetOutput()->Quantile(0,0.01);
  double maxModulus = modulusHistogram->GetOutput()->Quantile(0,0.98);

  ModulusFilterType::Pointer coherenceFilter = ModulusFilterType::New();
  coherenceFilter->SetInput(interferogram->GetOutput());

  typedef itk::ComplexToPhaseImageFilter<ImageType,ScalarImageType> PhaseFilterType;
  PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();
  phaseFilter->SetInput(interferogram->GetOutput());
  phaseFilter->Update();

  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef otb::Image<RGBPixelType, 2> RGBImageType;

  typedef AmplitudeCoherencyPhaseToRGBFunctor
      <ScalarPixelType,ScalarPixelType,ScalarPixelType,RGBPixelType> ColorMapFunctorType;
  typedef itk::TernaryFunctorImageFilter
      <ScalarImageType, ScalarImageType, ScalarImageType, RGBImageType, ColorMapFunctorType> ColorMapFilterType;
  ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
  colormapper->GetFunctor().SetMaximum(maxModulus);
  colormapper->GetFunctor().SetMinimum(minModulus);

  colormapper->SetInput1(modulusFilter->GetOutput());
  colormapper->SetInput2(coherenceFilter->GetOutput());
  colormapper->SetInput3(phaseFilter->GetOutput());

  typedef otb::StreamingImageFileWriter<RGBImageType> WriterRGBType;
  WriterRGBType::Pointer writerRGB = WriterRGBType::New();
  writerRGB->SetFileName(argv[4]);
  writerRGB->SetInput(colormapper->GetOutput());

  writerRGB->Update();

  return EXIT_SUCCESS;
}
