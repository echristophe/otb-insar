#ifndef __otbSARCoregistrationImageFilter_txx
#define __otbSARCoregistrationImageFilter_txx

#include "otbSARCoregistrationImageFilter.h"

#include "itkProgressReporter.h"
#include "otbWindowedSincInterpolateImageBlackmanFunction.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkExceptionObject.h"

namespace otb
{

/** SARCoregistrationImageFilter */

/**
 * Constructor
 */
template <class TInputImage, class TInterpolateFunction>
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::SARCoregistrationImageFilter()
 {
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfOutputs(1);

  // Default sizes
  m_TiePointsPerDim = 32;
  m_PatchSizePerDim = 128;
  m_SearchRadius.Fill(1);

  // Perform fine registration
  m_PerformFine = false;

  // Default sub-pixel precision
  m_SubPixelAccuracy = 0.125;

  // Flags
  m_UseSpacing = true;

  // Default interpolator
  m_Interpolator = otb::ComplexInterpolateImageFunction<
		TInputImage, 
		otb::Function::BlackmanWindowFunction< PixelType::value_type >, 
		BoundaryConditionType, double>::New();

  // Grid Step
  m_GridStep.Fill(1);

  m_Transform = NULL;

  m_UseDEM = false;
 }

template <class TInputImage, class TInterpolateFunction>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::SetMasterInput( const TInputImage * image )
 {
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage *>( image ));
 }

template <class TInputImage, class TInterpolateFunction>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::SetSlaveInput( const TInputImage * image)
 {
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TInputImage *>( image ));
 }

template <class TInputImage, class TInterpolateFunction>
const TInputImage *
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::GetMasterInput()
 {
  if (this->GetNumberOfInputs()<1)
    {
    return 0;
    }
  return static_cast<const TInputImage *>(this->itk::ProcessObject::GetInput(0));
 }

template <class TInputImage, class TInterpolateFunction>
const TInputImage *
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::GetSlaveInput()
 {
  if (this->GetNumberOfInputs()<2)
    {
    return 0;
    }
  return static_cast<const TInputImage *>(this->itk::ProcessObject::GetInput(1));
 }

template <class TInputImage, class TInterpolateFunction>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::GenerateOutputInformation()
 {
  // Call superclass implementation
  Superclass::GenerateOutputInformation();

  // Retrieve output pointers
  TInputImage * outputPtr = this->GetOutput();
  
  // Update size and spacing according to grid step
  InputImageRegionType largestRegion  = outputPtr->GetLargestPossibleRegion();
  SizeType outputSize       = largestRegion.GetSize();
  SpacingType outputSpacing = outputPtr->GetSpacing();

  for(unsigned int dim = 0; dim < TInputImage::ImageDimension; ++dim)
    {
    outputSize[dim] /= m_GridStep[dim];
    outputSpacing[dim] *= m_GridStep[dim];
    }

  // Set spacing
  outputPtr->SetSpacing(outputSpacing);

  // Set largest region size
  largestRegion.SetSize(outputSize);
  outputPtr->SetLargestPossibleRegion(largestRegion);
 }

template <class TInputImage, class TInterpolateFunction>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::GenerateInputRequestedRegion()
 {
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  TInputImage * masterPtr  = const_cast< TInputImage * >( this->GetMasterInput());
  TInputImage * slavePtr = const_cast< TInputImage * >( this->GetSlaveInput());

  TInputImage * outputPtr = this->GetOutput();

  if ( !masterPtr || !slavePtr || !outputPtr )
    {
    return;
    }

  // get a copy of the master requested region (should equal the output
  // requested region)
  InputImageRegionType masterRequestedRegion, slaveRequestedRegion;
  masterRequestedRegion = outputPtr->GetRequestedRegion();

  // Apply grid step
  SizeType masterRequestedSize = masterRequestedRegion.GetSize();
  IndexType masterRequestedIndex = masterRequestedRegion.GetIndex();

  for(unsigned int dim = 0; dim < TInputImage::ImageDimension; ++dim)
      {
      masterRequestedSize [dim] *= m_GridStep[dim];
      masterRequestedIndex[dim] *= m_GridStep[dim];
      }

  masterRequestedRegion.SetSize(masterRequestedSize);
  masterRequestedRegion.SetIndex(masterRequestedIndex);

  // pad the input requested region by the operator patchSize / 2
  masterRequestedRegion.PadByRadius( m_PatchSizePerDim / 2 );


  // get a copy of the slave requested region (should equal the output
  // requested region)
  InputImageRegionType searchMasterRequestedRegion = masterRequestedRegion;
  searchMasterRequestedRegion.PadByRadius(m_SearchRadius);


  // Find corners of the search window
   IndexType ulIndex = searchMasterRequestedRegion.GetIndex();

   IndexType lrIndex;
   for(unsigned int dim = 0; dim < TInputImage::ImageDimension; ++dim)
     {
     lrIndex[dim]= searchMasterRequestedRegion.GetIndex()[dim]
                 + searchMasterRequestedRegion.GetSize()[dim]-1;
     }

   // Transform back into slave region index space
   IndexType slaveIndex;

   // Find requested region
   SizeType slaveSize;

   slaveIndex = masterRequestedIndex;
   slaveSize = masterRequestedSize;

   slaveRequestedRegion.SetIndex(slaveIndex);
   slaveRequestedRegion.SetSize(slaveSize);

  // crop the fixed region at the fixed's largest possible region
  if ( masterRequestedRegion.Crop(masterPtr->GetLargestPossibleRegion()))
    {
    masterPtr->SetRequestedRegion( masterRequestedRegion );
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.
    // store what we tried to request (prior to trying to crop)
    masterPtr->SetRequestedRegion( masterRequestedRegion );

    // build an exception
    itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
    std::ostringstream msg;
    msg << this->GetNameOfClass()
                << "::GenerateInputRequestedRegion()";
    e.SetLocation(msg.str().c_str());
    e.SetDescription("Requested region is (at least partially) outside the largest possible region of image 1.");
    e.SetDataObject(masterPtr);
    throw e;
    }

  // crop the moving region at the moving's largest possible region
  if ( slaveRequestedRegion.Crop(slavePtr->GetLargestPossibleRegion()))
    {
    slavePtr->SetRequestedRegion( slaveRequestedRegion );
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region). This case might happen so we do not throw any exception but
    // request a null region instead
    slaveSize.Fill(0);
    slaveRequestedRegion.SetSize(slaveSize);
    slaveIndex.Fill(0);
    slaveRequestedRegion.SetIndex(slaveIndex);

    // store what we tried to request (prior to trying to crop)
    slavePtr->SetRequestedRegion(slaveRequestedRegion);
    }
  return;
 }

template <class TInputImage, class TInterpolateFunction>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction>
::GenerateData()
 {
  // Allocate outputs
  this->AllocateOutputs();

  // Get the image pointers
  const TInputImage * masterPtr = this->GetMasterInput();
  const TInputImage * slavePtr = this->GetSlaveInput();
  //TInputImage * outputPtr = this->GetOutput();
  
  // Wire currentMetric
  m_Interpolator->SetInputImage(this->GetSlaveInput());

  // Prepare Generic RS Transform
  m_Transform = TransformType::New();
  m_Transform->SetInputKeywordList(masterPtr->GetImageKeywordlist());
  m_Transform->SetOutputKeywordList(slavePtr->GetImageKeywordlist());
  if( m_UseDEM == true)
  {
	  m_Transform->SetDEMDirectory(m_DEMDir);
  }
  m_Transform->InstanciateTransform();

  // Define grid points from master image with requested number of Tie Points
  SizeType masterSize = masterPtr->GetRequestedRegion().GetSize();

  typename PointSetSourceType::Pointer pointSet = PointSetSourceType::New();
  typename PointSetType::PointType minPoint, maxPoint;

  SizeType patchSize;

  for(unsigned int i = 0; i < TInputImage::ImageDimension; i++)
  {
	  maxPoint[i] = masterSize[i] - ((masterSize[i] - 1) % m_TiePointsPerDim) - 1;
	  minPoint[i] = maxPoint[i] / m_TiePointsPerDim;
	  
	  // Set patch size for region	  
	  patchSize[i] = m_PatchSizePerDim;  
  }

  pointSet->SetMinPoint(minPoint);
  pointSet->SetMaxPoint(maxPoint);
  pointSet->SetNumberOfPoints(m_TiePointsPerDim);
  pointSet->Update();

  // Get reference to points from grid
  typename PointSetSourceType::PointsContainerPointer points;
  points = pointSet->GetOutput()->GetPoints();

  // Define estimator for affine transform
  typename EstimateFilterType::Pointer estimate = EstimateFilterType::New();

  // support progress methods/callbacks
  itk::ProgressReporter progress(this, 0, points->Size());

  // Walk the grid point by point
  typename PointsContainerType::ConstIterator gridIt = points->Begin();
  while (gridIt != points->End() )
    {
		PointType masterPoint = gridIt.Value();

		PointType slavePoint = m_Transform->TransformPoint(masterPoint);
		
		IndexType masterIndex;
		masterPtr->TransformPhysicalPointToIndex(masterPoint, masterIndex);
		
		TInputImage::PointType masterOrigin = masterPtr->GetOrigin();

		IndexType slaveIndex;
		slavePtr->TransformPhysicalPointToIndex(slavePoint, slaveIndex);

		TInputImage::PointType slaveOrigin = slavePtr->GetOrigin();

		IndexType currentMasterIndex;
		IndexType currentSlaveIndex;

		for(unsigned int i = 0; i < TInputImage::ImageDimension; i++)
		{
			currentMasterIndex[i] = masterIndex[i] - (m_PatchSizePerDim / 2) + (int)masterOrigin[i];
			currentSlaveIndex[i] = slaveIndex[i] - (m_PatchSizePerDim / 2) + (int)slaveOrigin[i];
		}

		RegionType masterCurrentRegion;
		masterCurrentRegion.SetIndex(currentMasterIndex);
		masterCurrentRegion.SetSize(patchSize);

		RegionType slaveCurrentRegion;
		slaveCurrentRegion.SetIndex(currentSlaveIndex);
		slaveCurrentRegion.SetSize(patchSize);

		if( !(slavePtr->GetLargestPossibleRegion().IsInside(slaveCurrentRegion)) || !(masterPtr->GetLargestPossibleRegion().IsInside(masterCurrentRegion)))
		{
			++gridIt;
			progress.CompletedPixel();

			continue;
		}

		// COARSE registration

		// Prepare extraction regions from images
		ExtractFilterType::Pointer masterExtract = ExtractFilterType::New();
		masterExtract->SetInput(masterPtr);
		masterExtract->SetExtractionRegion(masterCurrentRegion);

		ExtractFilterType::Pointer slaveExtract = ExtractFilterType::New();
		slaveExtract->SetInput(slavePtr);
		slaveExtract->SetExtractionRegion(slaveCurrentRegion);

		// Pad extraction regions
		SizeType padSize;
		padSize.Fill(m_PatchSizePerDim/2);

		PadFilterType::Pointer masterPad = PadFilterType::New();
		masterPad->SetInput(masterExtract->GetOutput());
		masterPad->SetPadBound(padSize);

		PadFilterType::Pointer slavePad = PadFilterType::New();
		slavePad->SetInput(slaveExtract->GetOutput());
		slavePad->SetPadBound(padSize);

		// Direct FFT on padded regions
		FFTType::Pointer fft = FFTType::New();
		fft->SetInput(masterPad->GetOutput());
		fft->Update();

		FFTOutputImageType::Pointer masterFFTImage = fft->GetOutput();

		fft = FFTType::New();
		fft->SetInput(slavePad->GetOutput());
		fft->Update();

		FFTOutputImageType::Pointer slaveFFTImage = fft->GetOutput();

		// Calculate complex conjugate product
		ConjugateProductFilterType::Pointer conjProduct = ConjugateProductFilterType::New();
		conjProduct->SetInput1(masterFFTImage);
		conjProduct->SetInput2(slaveFFTImage);

		// Pseudo normalize
		ModulusFilterType::Pointer conjProdModulus = ModulusFilterType::New();
		conjProdModulus->SetInput(conjProduct->GetOutput());

		DivideFilterType::Pointer conjProdNormalize = DivideFilterType::New();
		conjProdNormalize->SetInput1(conjProduct->GetOutput());
		conjProdNormalize->SetInput2(conjProdModulus->GetOutput());

		// Inverse FFT on normalized coefficient
		fft = FFTType::New();
		fft->SetInput(conjProdNormalize->GetOutput());
		fft->SetTransformDirection(FFTType::INVERSE);
		fft->Update();

		FFTOutputImageType::Pointer invFFT = fft->GetOutput();

		// Shift 
		ShiftFilterType::Pointer shift = ShiftFilterType::New();
		shift->SetInput(invFFT);

		// Get modulus (real valued coefficient)
		ModulusFilterType::Pointer modulus = ModulusFilterType::New();
		modulus->SetInput(shift->GetOutput());
		modulus->Update();

		// Get position of maximum correlation
		MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
		minMax->SetImage(modulus->GetOutput());
		minMax->ComputeMaximum();

		if( minMax->GetMaximum() > m_CorrelationThreshold)
		{
			IndexType maximumIndex = minMax->GetIndexOfMaximum();

			for(unsigned int i = 0; i < TInputImage::ImageDimension; i++)
			{
				slaveIndex[i] = slaveIndex[i] - maximumIndex[i] + (m_PatchSizePerDim / 2);
			}

			slavePtr->TransformIndexToPhysicalPoint(slaveIndex, slavePoint);

			estimate->AddTiePoints(masterPoint, slavePoint);
		}

   

		// FINE registration
		if( m_PerformFine == true)
		{
			// Is this really needed?

		}



    // Store the offset

    // Update iterator
    ++gridIt;

    // Update progress
    progress.CompletedPixel();
    }

	if( estimate->GetTiePointsContainer().size() != 0 )
	{
		estimate->Compute();

		NormalizeZeroFrequencyType::Pointer normalizeZeroFrequency = NormalizeZeroFrequencyType::New();
		normalizeZeroFrequency->SetImage(masterPtr);
		normalizeZeroFrequency->Compute();

		m_Interpolator->SetInputImage(slavePtr);
		m_Interpolator->SetRadius(3);
		m_Interpolator->SetNormalizeZeroFrequency(normalizeZeroFrequency->GetNormalizeZeroFrequency());

		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetTransform(estimate->GetAffineTransform());
		resample->SetInterpolator(m_Interpolator);
		resample->SetInput(slavePtr);
		resample->SetOutputSize(masterSize);
		resample->SetOutputOrigin(masterPtr->GetOrigin());
		resample->SetOutputSpacing(masterPtr->GetSpacing());
		resample->Update();

		this->GraftOutput(resample->GetOutput());
	}
 }
} // end namespace otb

#endif
