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
#ifndef __otbComplexInterpolateImageFunction_txx
#define __otbComplexInterpolateImageFunction_txx
#include "otbComplexInterpolateImageFunction.h"
#include "vnl/vnl_math.h"
#include "otbMath.h"

#include <complex>

namespace otb
{

/** Constructor */
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::ComplexInterpolateImageFunction()
{
  m_WindowSize = 1;
  this->SetRadius(1);
  m_OffsetTable = NULL;
  m_WeightOffsetTable = NULL;
  m_TablesHaveBeenGenerated = false;
  m_NormalizeWeight =  false;
  m_ZeroDoppler = 0.0;
}

/** Destructor */
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::~ComplexInterpolateImageFunction()
{
  this->ResetOffsetTable();
}

/** Delete every tables. */
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::ResetOffsetTable()
{
  // Clear the offset table
  if (m_OffsetTable != NULL)
    {
    delete[] m_OffsetTable;
    m_OffsetTable = NULL;
    }

  // Clear the weights tales
  if (m_WeightOffsetTable != NULL)
    {
    for (unsigned int i = 0; i < m_OffsetTableSize; ++i)
      {
      delete[] m_WeightOffsetTable[i];
      }
    delete[] m_WeightOffsetTable;
    m_WeightOffsetTable = NULL;
    }
}

template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::SetRadius(unsigned int rad)
{
  //m_Radius = rad;
  this->GetFunction().SetRadius(rad);
  m_WindowSize = rad << 1;
  this->Modified();
}

template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::Modified()
{
  Superclass::Modified();
  m_TablesHaveBeenGenerated = false;

}

/** Initialize used tables*/
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::InitializeTables()
{
  // Compute the offset table size
  m_OffsetTableSize = 1;
  for (unsigned dim = 0; dim < ImageDimension; ++dim)
    {
    m_OffsetTableSize *= m_WindowSize;
    }

  // Allocate the offset table
  m_OffsetTable = new unsigned int[m_OffsetTableSize];

  // Allocate the weights tables
  m_WeightOffsetTable = new unsigned int *[m_OffsetTableSize];
  for (unsigned int i = 0; i < m_OffsetTableSize; ++i)
    {
    m_WeightOffsetTable[i] = new unsigned int[ImageDimension];
    }
}

/** Fill the weight offset table*/
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::FillWeightOffsetTable()
{
  // Initialize the neighborhood
  SizeType radius;
  radius.Fill(this->GetRadius());
  if (this->GetInputImage() != NULL)
    {
    IteratorType it = IteratorType(radius,  this->GetInputImage(), this->GetInputImage()->GetBufferedRegion());
    // Compute the offset tables (we ignore all the zero indices
    // in the neighborhood)
    unsigned int iOffset = 0;
    int          empty = static_cast<int>(this->GetRadius());

    for (unsigned int iPos = 0; iPos < it.Size(); ++iPos)
      {
      // Get the offset (index)
      typename IteratorType::OffsetType off = it.GetOffset(iPos);

      // Check if the offset has zero weights
      bool nonzero = true;
      for (unsigned int dim = 0; dim < ImageDimension; ++dim)
        {
        if (off[dim] == -empty)
          {
          nonzero = false;
          break;
          }
        }
      // Only use offsets with non-zero indices
      if (nonzero)
        {
        // Set the offset index
        m_OffsetTable[iOffset] = iPos;

        // Set the weight table indices
        for (unsigned int dim = 0; dim < ImageDimension; ++dim)
          {
          m_WeightOffsetTable[iOffset][dim] = off[dim] + this->GetRadius() - 1;
          }
        // Increment the index
        iOffset++;
        }
      }
    }
  else
    {
    itkExceptionMacro(<< "An input has to be set");
    }
}

/** Initialize tables: need to be call explicitely */
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::Initialize()
{
  // Delete existing tables
  this->ResetOffsetTable();
  // Tables initialization
  this->InitializeTables();
  // fill the weigth table
  this->FillWeightOffsetTable();
  m_TablesHaveBeenGenerated = true;
}

/** Evaluate at image index position */
template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
typename ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>::OutputType
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::EvaluateAtContinuousIndex(const ContinuousIndexType& index) const
{
  IndexType baseIndex;

  for (unsigned int dim = 0; dim < ImageDimension; ++dim)
	{
	baseIndex[dim] = std::floor(index[dim]);
	}

  // Position the neighborhood at the index of interest
  SizeType radius;
  radius.Fill(this->GetRadius());
  IteratorType nit = IteratorType(radius, this->GetInputImage(), this->GetInputImage()->GetBufferedRegion());
  nit.SetLocation(baseIndex);

  RealType sumFunction = 0.;
  RealType resultValue = 0.;

  for(unsigned int elt = 0 ; elt < nit.Size(); ++elt) 
	{
		IteratorType::OffsetType offset = nit.GetOffset(elt);
		RealType  valueFunction = 1.0;
		RealType PhaseShift(1.0,0.0);
		for(unsigned int dim = 0; dim < ImageDimension; ++dim)
			{
				ScalarRealType delta = otb::CONST_2PI * m_ZeroDoppler * offset[dim];
				RealType localPhaseShift(cos(delta),sin(delta));
				PhaseShift *= localPhaseShift;
				RealType valueTmp = m_Function(offset[dim]); 
				valueFunction *=  valueTmp;
			}

		sumFunction = sumFunction + valueFunction;
		RealType pixelValue = static_cast<RealType>(nit.GetPixel(elt));
		resultValue += valueFunction * pixelValue * PhaseShift; 
  }


  // Return the interpolated value
  //OutputType result;
  //result = this->ConvertValue(resultValue);
  //return result;
  return 0.0;
}

template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
typename ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>::OutputType
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::ConvertValue(std::complex<ScalarRealType> value)
{
	return static_cast<OutputType>(std::abs(value)); 
}

template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
typename ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>::OutputType
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::ConvertValue(ScalarRealType value)
{
	return static_cast<OutputType>(value); 
}

template<class TInputImage, class TFunction, class TBoundaryCondition, class TCoordRep>
void
ComplexInterpolateImageFunction<TInputImage, TFunction, TBoundaryCondition, TCoordRep>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} //namespace otb

#endif
