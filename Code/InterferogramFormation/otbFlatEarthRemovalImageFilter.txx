#ifndef __otbFlatEarthRemovalImageFilter_txx
#define __otbFlatEarthRemovalImageFilter_txx

#include "otbFlatEarthRemovalImageFilter.h"

#include "itkProgressReporter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkExceptionObject.h"

namespace otb
{

/** FlatEarthRemovalImageFilter */

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
FlatEarthRemovalImageFilter<TInputImage, TOutputImage>
::FlatEarthRemovalImageFilter()
 {
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfOutputs(1);

  // Default sizes
  m_PatchSizePerDim = 8;
  m_PadSizePerDim = 8;
 }

template <class TInputImage, class TOutputImage>
void
FlatEarthRemovalImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
 {
  // Call superclass implementation
  Superclass::GenerateOutputInformation();

  // Retrieve output pointers
  TOutputImage * outputPtr = this->GetOutput();
  
  // Update size and spacing according to grid step
  InputImageRegionType largestRegion  = outputPtr->GetLargestPossibleRegion();
  outputPtr->SetLargestPossibleRegion(largestRegion);
 }

template <class TInputImage, class TOutputImage>
void
FlatEarthRemovalImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
 {
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  TInputImage * inputPtr  = const_cast< TInputImage * >( this->GetInput());
  
  TOutputImage * outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // get a copy of the master requested region (should equal the output
  // requested region)
  InputImageRegionType inputRequestedRegion;
  inputRequestedRegion = outputPtr->GetRequestedRegion();

  // pad the input requested region by the operator patchSize / 2
  inputRequestedRegion.PadByRadius( m_PatchSizePerDim / 2 );

  // crop the fixed region at the fixed's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
    {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.
    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );

    // build an exception
    itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
    std::ostringstream msg;
    msg << this->GetNameOfClass()
                << "::GenerateInputRequestedRegion()";
    e.SetLocation(msg.str().c_str());
    e.SetDescription("Requested region is (at least partially) outside the largest possible region of image 1.");
    e.SetDataObject(inputPtr);
    throw e;
    }
  
  return;
 }

template <class TInputImage, class TOutputImage>
void
FlatEarthRemovalImageFilter<TInputImage, TOutputImage>
::GenerateData()
 {
  // Allocate outputs
  this->AllocateOutputs();

  // Get the image pointers
  const TInputImage * inputPtr = this->GetInput();
  TOutputImage * outputPtr = this->GetOutput();
  
  SizeType inputSize = inputPtr->GetRequestedRegion().GetSize();

  // support progress methods/callbacks
  //itk::ProgressReporter progress(this, 0, in);

  SubsampledRegionIteratorType inIt(inputPtr, inputPtr->GetRequestedRegion());

  // Walk the image
  IndexType offset;
  offset[0] = m_PatchSizePerDim;
  offset[1] = m_PatchSizePerDim;
  inIt.SetSubsampleFactor(offset);
  inIt.GoToBegin();
  while (!inIt.IsAtEnd() )
    {
		IndexType currentIndex = inIt.GetIndex();

		SizeType patchSize;
		for(unsigned int i = 0; i < TInputImage::ImageDimension; i++)
		{
			patchSize[i] = m_PatchSizePerDim;
		}

		RegionType currentRegion;
		currentRegion.SetIndex(currentIndex);
		currentRegion.SetSize(patchSize);

		currentRegion.Crop(inputPtr->GetRequestedRegion());
		
		// Prepare extraction region from image
		typename ExtractFilterType::Pointer extract = ExtractFilterType::New();
		extract->SetInput(inputPtr);
		extract->SetExtractionRegion(currentRegion);
		extract->Update();

		typename TInputImage::Pointer extractImage = extract->GetOutput();

		typedef typename otb::ImageFileWriter<TInputImage> InputWriterType;
		typename InputWriterType::Pointer extractWriter = InputWriterType::New();
		extractWriter->SetInput(extract->GetOutput());
		extractWriter->SetFileName("flattests/RS2extract.tif");
		//extractWriter->Update();

		typename TInputImage::SizeType paddsize;
        paddsize.Fill(m_PatchSizePerDim/2);

		typename PadFilterType::Pointer pad = PadFilterType::New();
		pad->SetInput(extract->GetOutput());
		pad->SetPadBound(paddsize);
		pad->SetConstant(0.0);

		// Direct FFT on padded regions
		typename FFTType::Pointer fft = FFTType::New();
		//fft->SetInput(extract->GetOutput());
		fft->SetInput(pad->GetOutput());
		fft->Update();

		typename FFTOutputImageType::Pointer fftImage = fft->GetOutput();

		// Get modulus
		typename ModulusFilterType::Pointer modulus = ModulusFilterType::New();
		modulus->SetInput(fftImage);
		modulus->Update();

		// Get position of peak
		typename MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
		minMax->SetImage(modulus->GetOutput());
		minMax->ComputeMaximum();

		IndexType maximumIndex = minMax->GetIndexOfMaximum();

		// Get original phase image for block
		typename PhaseFilterType::Pointer phase = PhaseFilterType::New();
		phase->SetInput(extract->GetOutput());
		phase->Update();

		typename TOutputImage::Pointer originalPhaseImage = phase->GetOutput();

		double preciseIndex[2];

		for(unsigned int i = 0; i < TInputImage::ImageDimension; i++)
		{
			preciseIndex[i] = ((double)(maximumIndex[i] + m_PadSizePerDim/2) / 2);			
			maximumIndex[i] = (typename IndexType::IndexValueType)((maximumIndex[i] + m_PadSizePerDim/2) / 2);			
		}

		typename TOutputImage::PixelType phasePeakValue = originalPhaseImage->GetPixel(maximumIndex);

		typename FlatEarthPhaseCalculationType::Pointer flatEarthPhaseCalculate = FlatEarthPhaseCalculationType::New();
		flatEarthPhaseCalculate->SetInput(originalPhaseImage);
		flatEarthPhaseCalculate->GetFunctor().SetFringePhase(phasePeakValue);
		flatEarthPhaseCalculate->GetFunctor().SetRangeRate(preciseIndex[0]);
		flatEarthPhaseCalculate->GetFunctor().SetAzimuthRate(preciseIndex[1]);
		flatEarthPhaseCalculate->Update();

		typename TInputImage::Pointer estimatedPhaseImage = flatEarthPhaseCalculate->GetOutput();

		typename MultiplyFilterType::Pointer multiply = MultiplyFilterType::New();
		multiply->SetInput1(extract->GetOutput());
		multiply->SetInput2(estimatedPhaseImage);
		multiply->Update();

		typename TInputImage::Pointer flatRemovedImage = multiply->GetOutput();
		OutRegionIteratorType outIt(outputPtr, currentRegion);
		InRegionIteratorType flIt(flatRemovedImage, flatRemovedImage->GetRequestedRegion());

		for( outIt.GoToBegin(), flIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt, ++flIt )
		{
			PixelType cplxValue = flIt.Get();
			outIt.Set(atan2(cplxValue.imag(), cplxValue.real()));
		}

    // Update iterator
    ++inIt;

    // Update progress
    //progress.CompletedPixel();
    }
 }
} // end namespace otb

#endif
