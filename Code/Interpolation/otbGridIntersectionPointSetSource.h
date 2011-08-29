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
