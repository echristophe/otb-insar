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
