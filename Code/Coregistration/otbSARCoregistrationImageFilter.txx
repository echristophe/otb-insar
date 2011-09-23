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
#ifndef __otbSARCoregistrationImageFilter_txx
#define __otbSARCoregistrationImageFilter_txx

#include "otbSARCoregistrationImageFilter.h"

#include "itkProgressReporter.h"
#include "otbWindowedSincInterpolateImageBlackmanFunction.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkExceptionObject.h"

namespace otb
{
/**
 * Constructor
 */
template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, T0utputCorrelation, TOutputDeformationField>
::SARCoregistrationImageFilter()
 {
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfOutputs(2);
  this->SetNthOutput(1, TOutputDeformationField::New());

  // Default sizes
  m_TiePointsPerDim = 32;
  m_PatchSizePerDim = 128;
  m_SearchRadius.Fill(1);

  // Default sub-pixel precision
  m_SubPixelAccuracy = 0.1;

  // Flags
  m_UseSpacing = true;

  // Default interpolator
  m_Interpolator = otb::ComplexInterpolateImageFunction<TInputImage, otb::Function::BlackmanWindowFunction< PixelType::value_type >, BoundaryConditionType, double>::New();

  // Grid Step
  m_GridStep.Fill(1);

  m_Transform = NULL;
 }

template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, T0utputCorrelation, TOutputDeformationField>
::SetMasterInput( const TInputImage, TInterpolateFunction * image )
 {
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage, TInterpolateFunction *>( image ));
 }

template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, T0utputCorrelation, TOutputDeformationField>
::SetSlaveInput( const TInputImage, TInterpolateFunction * image)
 {
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TInputImage, TInterpolateFunction *>( image ));
 }

template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
const TInputImage, TInterpolateFunction *
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, T0utputCorrelation, TOutputDeformationField>
::GetMasterInput()
 {
  if (this->GetNumberOfInputs()<1)
    {
    return 0;
    }
  return static_cast<const TInputImage, TInterpolateFunction *>(this->itk::ProcessObject::GetInput(0));
 }

template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
const TInputImage, TInterpolateFunction *
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, T0utputCorrelation, TOutputDeformationField>
::GetSlaveInput()
 {
  if (this->GetNumberOfInputs()<2)
    {
    return 0;
    }
  return static_cast<const TInputImage, TInterpolateFunction *>(this->itk::ProcessObject::GetInput(1));
 }

template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
TOutputDeformationField *
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, T0utputCorrelation, TOutputDeformationField>
::GetOutputDeformationField()
 {
  if (this->GetNumberOfOutputs()<2)
    {
    return 0;
    }
  return static_cast<TOutputDeformationField *>(this->itk::ProcessObject::GetOutput(1));
 }

template <class TInputImage, class TInterpolateFunction, class TOutputCorrelation, class TOutputDeformationField>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, TOutputCorrelation, TOutputDeformationField>
::GenerateOutputInformation()
 {
  // Call superclass implementation
  Superclass::GenerateOutputInformation();

  // Retrieve output pointers
  TOutputCorrelation * outputPtr = this->GetOutput();
  TOutputDeformationField *outputFieldPtr = this->GetOutputDeformationField();

  // Update size and spacing according to grid step
  InputImageRegionType largestRegion  = outputPtr->GetLargestPossibleRegion();
  SizeType outputSize       = largestRegion.GetSize();
  SpacingType outputSpacing = outputPtr->GetSpacing();

  for(unsigned int dim = 0; dim < TOutputCorrelation::ImageDimension; ++dim)
    {
    outputSize[dim] /= m_GridStep[dim];
    outputSpacing[dim] *= m_GridStep[dim];
    }

  // Set spacing
  outputPtr->SetSpacing(outputSpacing);
  outputFieldPtr->SetSpacing(outputSpacing);

  // Set largest region size
  largestRegion.SetSize(outputSize);
  outputPtr->SetLargestPossibleRegion(largestRegion);
  outputFieldPtr->SetLargestPossibleRegion(largestRegion);
 }

template <class TInputImage, class TInterpolateFunction, class TOutputCorrelation, class TOutputDeformationField>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, TOutputCorrelation, TOutputDeformationField>
::GenerateInputRequestedRegion()
 {
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  TInputImage, TInterpolateFunction * masterPtr  = const_cast< TInputImage, TInterpolateFunction * >( this->GetFixedInput());
  TInputImage, TInterpolateFunction * slavePtr = const_cast< TInputImage, TInterpolateFunction * >( this->GetMovingInput());

  TOutputCorrelation * outputPtr = this->GetOutput();

  if ( !masterPtr || !slavePtr || !outputPtr )
    {
    return;
    }

  // get a copy of the fixed requested region (should equal the output
  // requested region)
  InputImageRegionType masterRequestedRegion, slaveRequestedRegion;
  masterRequestedRegion = outputPtr->GetRequestedRegion();

  // Apply grid step
  SizeType masterRequestedSize = masterRequestedRegion.GetSize();
  IndexType masterRequestedIndex = masterRequestedRegion.GetIndex();

  for(unsigned int dim = 0; dim < TOutputCorrelation::ImageDimension; ++dim)
      {
      masterRequestedSize [dim] *= m_GridStep[dim];
      masterRequestedIndex[dim] *= m_GridStep[dim];
      }

  masterRequestedRegion.SetSize(masterRequestedSize);
  masterRequestedRegion.SetIndex(masterRequestedIndex);

  // pad the input requested region by the operator patchSize / 2
  masterRequestedRegion.PadByRadius( m_PatchSizePerDim / 2 );


  // get a copy of the moving requested region (should equal the output
  // requested region)
  InputImageRegionType searchMasterRequestedRegion = masterRequestedRegion;
  searchMasterRequestedRegion.PadByRadius(m_SearchRadius);


  // Find corners of the search window
   IndexType ulIndex = searchMasterRequestedRegion.GetIndex();

   IndexType lrIndex;
   for(unsigned int dim = 0; dim < TInputImage, TInterpolateFunction::ImageDimension; ++dim)
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

template <class TInputImage, class TInterpolateFunction, class TOutputCorrelation, class TOutputDeformationField>
void
SARCoregistrationImageFilter<TInputImage, TInterpolateFunction, TOutputCorrelation, TOutputDeformationField>
::GenerateData()
 {
  // Allocate outputs
  this->AllocateOutputs();

  // Get the image pointers
  const TInputImage, TInterpolateFunction * masterPtr = this->GetMasterInput();
  const TInputImage, TInterpolateFunction * slavePtr = this->GetSlaveInput();
  TOutputCorrelation * outputPtr = this->GetOutput();
  TOutputDeformationField * outputDfPtr = this->GetOutputDeformationField();

  // Wire currentMetric
  m_Interpolator->SeInputImage(this->GetSlaveInput());
  
  /** Output iterators */
  itk::ImageRegionIteratorWithIndex<TOutputCorrelation> outputIt(outputPtr, outputPtr->GetRequestedRegion());
  itk::ImageRegionIterator<TOutputDeformationField> outputDfIt(outputDfPtr, outputPtr->GetRequestedRegion());
  outputIt.GoToBegin();
  outputDfIt.GoToBegin();

  // support progress methods/callbacks
  itk::ProgressReporter progress(this, 0, outputPtr->GetRequestedRegion().GetNumberOfPixels());

  // GoToBegin
  outputIt.GoToBegin();
  outputDfIt.GoToBegin();

  // Correl, max correl, maxPosition
  double currentMetric, optMetric;

  // Optimal translation parameters
  typename TranslationType::ParametersType params(2), optParams(2), tmpOptParams(2);

  // Final deformation value
  DeformationValueType deformationValue;
  deformationValue[0] = m_InitialOffset[0];
  deformationValue[1] = m_InitialOffset[1];

  // Local initial offset: enable the possibility of a different initial offset for each pixel
  SpacingType localOffset = m_InitialOffset;

  // Get fixed image spacing
  SpacingType fixedSpacing = fixedPtr->GetSpacing();

  // Walk the images
  while (!outputIt.IsAtEnd() && !outputDfIt.IsAtEnd() )
    {
    
    
    // Apply grid step
    IndexType currentIndex = outputIt.GetIndex();
    for(unsigned int dim = 0; dim < TInputImage, TInterpolateFunction::ImageDimension; ++dim)
      {
      currentIndex[dim] *= m_GridStep[dim];
      }
    

    // Compute the local offset if required (and the transform was specified)
    if (m_Transform.IsNotNull())
      {
      PointType inputPoint, outputPoint;
      for(unsigned int dim = 0; dim < TInputImage, TInterpolateFunction::ImageDimension; ++dim)
        {
        inputPoint[dim] = currentIndex[dim];
        }
      outputPoint = m_Transform->TransformPoint(inputPoint);
      for(unsigned int dim = 0; dim < TInputImage, TInterpolateFunction::ImageDimension; ++dim)
        {
        localOffset[dim] = outputPoint[dim] - inputPoint[dim]; //FIXME check the direction
        }
      }

   

    // Dichotomic sub-pixel
    SpacingType subPixelSpacing = fixedSpacing;
    while(subPixelSpacing[0] > m_SubPixelAccuracy || subPixelSpacing[1] > m_SubPixelAccuracy)
      {
      // Perform 1 step of dichotomic search
      subPixelSpacing /= 2.;

      // Store last opt params
      tmpOptParams = optParams;

      for(int i = -1; i <= 1; i+=2)
        {
        for(int j = -1; j <= 1; j+=2)
          {
          params = tmpOptParams;
          params[0] += static_cast<double>(i*subPixelSpacing[0]);
          params[1] += static_cast<double>(j*subPixelSpacing[1]);

          try
          {
            // compute currentMetric
            currentMetric = m_Metric->GetValue(params);

            // Check for maximum
            if((m_Minimize && (currentMetric < optMetric)) || (!m_Minimize && (currentMetric > optMetric)))
              {
              optMetric = currentMetric;
              optParams = params;
              }
          }
          catch(itk::ExceptionObject& err)
          {
            itkWarningMacro(<<err.GetDescription());

          }
          }
        }
      }

    // Store the offset and the correlation value
    outputIt.Set(optMetric);
    if(m_UseSpacing)
      {
      deformationValue[0] = optParams[0];
      deformationValue[1] = optParams[1];
      }
    else
      {
      deformationValue[0] = optParams[0]/fixedSpacing[0];
      deformationValue[1] = optParams[1]/fixedSpacing[1];
      }
    outputDfIt.Set(deformationValue);
    // Update iterators
    ++outputIt;
    ++outputDfIt;

    // Update progress
    progress.CompletedPixel();
    }
 }
} // end namespace otb

#endif
