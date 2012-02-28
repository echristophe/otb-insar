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
#ifndef __otbParallelPerpendicularBaselineFunctor_cxx
#define __otbParallelPerpendicularBaselineFunctor_cxx

#include "otbParallelPerpendicularBaselineFunctor.h"
#include <vnl/vnl_cross.h>

namespace otb
{ 

namespace Functor {


ParallelPerpendicularBaselineFunctor::OutputType 
ParallelPerpendicularBaselineFunctor
::GetParallelBaseline() const
{
	
	VectorType parallelVector;
	parallelVector = m_RefPoint;
	parallelVector(1) = 0.0;

	return dot_product(this->GetRTNBaseline(),parallelVector.normalize());
}

ParallelPerpendicularBaselineFunctor::OutputType 
ParallelPerpendicularBaselineFunctor
::GetPerpendicularBaseline() const
{

	VectorType parallelVector;
	parallelVector = m_RefPoint;
	parallelVector(1) = 0.0;

	VectorType tangentialVector;
	tangentialVector.fill(0.0);
	tangentialVector(1) = 1.0;

	VectorType perpendicularVector;
	perpendicularVector = vnl_cross_3d(parallelVector,tangentialVector);

	return dot_product(this->GetRTNBaseline(),perpendicularVector.normalize());
	itkExceptionMacro("GetPerpendicularBaseline() not implemented");
}


ParallelPerpendicularBaselineFunctor::OutputType 
ParallelPerpendicularBaselineFunctor
::GetBaseline(BaselineCalculusEnumType map)
{
  switch( map )
    {
    case Parallel:
      {
      return this->GetParallelBaseline();
      break;
      }
    case Perpendicular:
      {
      return this->GetPerpendicularBaseline();
      break;
      }
    }
	return 0;
}


void
ParallelPerpendicularBaselineFunctor
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "RefPoint : " << m_RefPoint << std::endl;
}


} // end namespace functor
} // end namespace otb

#endif
