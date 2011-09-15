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
#ifndef __otbImageNormalizeZeroFrequencyCalculator_txx
#define __otbImageNormalizeZeroFrequencyCalculator_txx
#include "otbImageNormalizeZeroFrequencyCalculator.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkFFTComplexToComplexImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkBinaryThresholdImageFilter.h"

namespace otb
{ 


class ITK_EXPORT InvalidImageMomentsError : public itk::ExceptionObject
{
public:
  /**
   * Constructor. Needed to ensure the exception object can be copied.
   */
  InvalidImageMomentsError(const char *file, unsigned int lineNumber) : ExceptionObject(file, lineNumber) { this->SetDescription("No valid image moments are availble.");}

  /**
   * Constructor. Needed to ensure the exception object can be copied.
   */
  InvalidImageMomentsError(const std::string& file, unsigned int lineNumber) : ExceptionObject(file, lineNumber) { this->SetDescription("No valid image moments are availble.");}  
  
  itkTypeMacro(InvalidImageMomentsError, ExceptionObject);
};

  
//----------------------------------------------------------------------
// Construct without computing moments
template<class TImage>
ImageNormalizeZeroFrequencyCalculator<TImage>::ImageNormalizeZeroFrequencyCalculator(void) 
{
  m_Valid = false;
  m_Image = NULL;
  m_NormalizeZeroFrequency.Fill(itk::NumericTraits<ITK_TYPENAME VectorType::ValueType>::Zero);
}

//----------------------------------------------------------------------
// Destructor
template<class TImage>
ImageNormalizeZeroFrequencyCalculator<TImage>::
~ImageNormalizeZeroFrequencyCalculator()
{
}

template<class TInputImage>
void
ImageNormalizeZeroFrequencyCalculator<TInputImage>
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
  os << indent << "Valid: " << m_Valid << std::endl;
  os << indent << "Normalize Zero Frequency : " << m_NormalizeZeroFrequency << std::endl;
}

//----------------------------------------------------------------------
// Compute Normalize Zero Frequency for a new or modified image
template<class TImage>
void
ImageNormalizeZeroFrequencyCalculator<TImage>::
Compute()
{
  m_NormalizeZeroFrequency = itk::NumericTraits<ScalarType>::Zero;
  typedef typename ImageType::IndexType IndexType;

  if( !m_Image ) 
    {
    return;
    }

  typedef itk::FFTComplexToComplexImageFilter<ValueType,ImageDimension> FFTImageType;

  typedef itk::Image<typename FFTImageType::OutputImagePixelType,ImageDimension> FFTOutputImageType;
  typedef itk::ComplexToModulusImageFilter<FFTOutputImageType,ScalarImageType>   ComplexToModulusImageType;
  typedef itk::MinimumMaximumImageCalculator<ScalarImageType> MaxImageCalulatorType;
  typedef itk::BinaryThresholdImageFilter< ScalarImageType, ScalarImageType > ThresholdImageFilterType;

  typename FFTImageType::Pointer fftImage = FFTImageType::New(); 
  typename ComplexToModulusImageType::Pointer modulusImage = ComplexToModulusImageType::New(); 

  fftImage->SetInput(m_Image.GetPointer());

  fftImage->SetTransformDirection(FFTImageType::DIRECT);
  typename FFTOutputImageType::SpacingType spacingFFT;

  fftImage->GetOutput()->UpdateOutputInformation();

  for(unsigned int dim = 0 ; dim < ImageDimension ; ++dim)
	{	
	spacingFFT[dim] = 1.0 / fftImage->GetOutput()->GetLargestPossibleRegion ().GetSize()[dim];
	}	

  modulusImage->SetInput(fftImage->GetOutput());
  modulusImage->Update();

  typename MaxImageCalulatorType::Pointer  maxCalculator = MaxImageCalulatorType::New();

  maxCalculator->SetImage(modulusImage->GetOutput());
  maxCalculator->ComputeMaximum();
  	  
  ScalarPixelType maximumValue;

  maximumValue = maxCalculator->GetMaximum();
  //std::cout << "Maximum : "  << maximumValue << std::endl;
  //std::cout << "IndexMax : " << maxCalculator->GetIndexOfMaximum() << std::endl;

  typename ThresholdImageFilterType::Pointer thresholdImage = ThresholdImageFilterType::New();
  thresholdImage->SetInput(modulusImage->GetOutput());
  thresholdImage->SetOutsideValue(0.0);
  thresholdImage->SetInsideValue(1.0);
  thresholdImage->SetLowerThreshold(maximumValue/2.0 );
  thresholdImage->SetUpperThreshold(maximumValue );
  thresholdImage->Update();

  itk::ImageRegionConstIteratorWithIndex< ScalarImageType > it(	thresholdImage->GetOutput(),
																thresholdImage->GetOutput()->GetRequestedRegion() ); 
  
  IndexType averageIndex;
  averageIndex.Fill(0);
  double nbIndex =0.0;

  it.Begin();

  while( !it.IsAtEnd() )
    {
    IndexType indexPosition = it.GetIndex();

	double value =it.Get();
	if(value >0.0)
		{
		for(unsigned int dim = 0 ; dim < ImageDimension ; ++dim)
			{	
			averageIndex[dim] += (indexPosition[dim]) ; 	
			}	
		++nbIndex;
		}
    ++it;
	}
  //std::cout << "averageIndex : " << averageIndex[0] << " , "<< averageIndex[1] << std::endl;


 for(unsigned int dim = 0 ; dim < ImageDimension ; ++dim)
	{	
	m_NormalizeZeroFrequency[dim] = static_cast<ScalarType>(averageIndex[dim] / nbIndex * spacingFFT[dim]); 	
	}	


  /* Remember that the moments are valid */
  m_Valid = 1;

}


//---------------------------------------------------------------------
// Get Normalize Zero Frequency
template<class TImage>
typename ImageNormalizeZeroFrequencyCalculator<TImage>::VectorType
ImageNormalizeZeroFrequencyCalculator<TImage>::
GetNormalizeZeroFrequency() const
{
  if (!m_Valid) 
    {
    itkExceptionMacro( << "GetNormalizeZeroFrequency() invoked, but the normalize Zero frequency have not been computed. Call Compute() first.");
    }
  return m_NormalizeZeroFrequency;
}

} // end namespace otb

#endif
