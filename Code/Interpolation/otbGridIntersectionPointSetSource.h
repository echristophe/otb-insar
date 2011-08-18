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

#ifndef __otbGridIntersectionPointSetSource_h
#define __otbGridIntersectionPointSetSource_h

#include "otbPointSetSource.h"
#include "itkVariableLengthVector.h"

namespace otb
{

/** \class GridIntersectionPointSetSource
 *  \brief This class generate a regular grid intersetion point set
 *
 * \ingroup DataSources
 */

template <class TOutputPointSet>
class ITK_EXPORT GridIntersectionPointSetSource
  : public PointSetSource<TOutputPointSet>
{

public:
  /** Standard class typedefs. */
  typedef GridIntersectionPointSetSource  Self;
  typedef PointSetSource<TOutputPointSet> Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods). */
  itkTypeMacro(GridIntersectionPointSetSource, PointSetSource)

  /** Some convenient typedefs. */
  typedef itk::DataObject::Pointer                     DataObjectPointer;
  typedef TOutputPointSet                              OutputPointSetType;
  typedef typename OutputPointSetType::Pointer         OutputPointSetPointer;
  typedef typename OutputPointSetType::PointsContainer PointsContainerType;
  typedef typename PointsContainerType::Pointer        PointsContainerPointer;
  typedef typename OutputPointSetType::PointType       PointType;

  /** VariableLengthVectorType typedef support*/
  typedef itk::VariableLengthVector<unsigned int> VariableLengthVectorType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(PointDimension, unsigned int, TOutputPointSet::PointDimension);

  itkSetMacro(MinPoint, PointType)
  itkGetMacro(MinPoint, PointType)
  itkSetMacro(MaxPoint, PointType)
  itkGetMacro(MaxPoint, PointType)

  void SetNumberOfPoints(unsigned int numberOfPoints)
  {
	m_NumberOfPoints.Fill(numberOfPoints);
  }

  void SetNumberOfPoints(VariableLengthVectorType value)
  {
	  if(value.Size() != m_NumberOfPoints.Size())
	  {
		itkExceptionMacro(<< "The agrument dimension of the SetNumberOfPoints() method must be the same dimension as the Point");
	  }
	  m_NumberOfPoints = value ;
  }



protected:
  GridIntersectionPointSetSource();
  virtual ~GridIntersectionPointSetSource() {}

  virtual void GenerateData(void);

private:
  GridIntersectionPointSetSource(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  VariableLengthVectorType m_NumberOfPoints;

  PointType m_MinPoint;
  PointType m_MaxPoint;

};

}

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbGridIntersectionPointSetSource.txx"
#endif

#endif
