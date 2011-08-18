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

#ifndef __otbGridIntersectionPointSetSource_txx
#define __otbGridIntersectionPointSetSource_txx

#include "otbGridIntersectionPointSetSource.h"

namespace otb
{

template<class TOutputPointSet>
GridIntersectionPointSetSource<TOutputPointSet>
::GridIntersectionPointSetSource()
{
  m_MinPoint.Fill(0.0);
  m_MaxPoint.Fill(1.0);
  m_NumberOfPoints.SetSize(PointDimension);
  m_NumberOfPoints.Fill(1);
}


template<class TOutputPointSet>
void
GridIntersectionPointSetSource<TOutputPointSet>
::GenerateData(void)
{
  OutputPointSetPointer outputPointSet = this->GetOutput();
  outputPointSet->Initialize();

  PointsContainerPointer outPoints = outputPointSet->GetPoints();

  unsigned int numberOfPoints = m_NumberOfPoints[0];

  for(unsigned int dim = 1 ; dim < PointDimension; ++dim)
	{
		numberOfPoints *= m_NumberOfPoints[dim];
	}

  outPoints->Reserve(numberOfPoints);

  typename PointsContainerType::Iterator outputPoint = outPoints->Begin();
  
  VariableLengthVectorType index;
  index.SetSize(PointDimension+1);
  index.Fill(1);
  
  PointType pointValue;
  for(unsigned int i = 0 ; i < m_NumberOfPoints[0]; ++i)
  {
	if(m_NumberOfPoints[0] > 1)
		{
			pointValue[0] = m_MinPoint[0] + (m_MaxPoint[0]-m_MinPoint[0])/(m_NumberOfPoints[0]-1) * i ;
		}
	  else
		{
			pointValue[0] = (m_MinPoint[0]+m_MaxPoint[0])/2 ;
		}
	for(unsigned int j = 0 ; j < m_NumberOfPoints[1]; ++j)
		{
		  if(m_NumberOfPoints[1] > 1)
			{
				pointValue[1] = m_MinPoint[1] + (m_MaxPoint[1]-m_MinPoint[1])/(m_NumberOfPoints[1]-1) * j ;
			}
		  else
			{
			pointValue[1] = (m_MinPoint[1]+m_MaxPoint[1])/2 ;
			}
		outputPoint.Value() = pointValue;
		++outputPoint;
		}
  }

}
}
#endif
