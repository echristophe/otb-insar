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
	
	VectorType parallelVector(3);
	parallelVector(0) = cos(m_LookDirection);
	parallelVector(1) = 0.0;
	parallelVector(2) = sin(m_LookDirection);
	
	return element_product(this->GetRTNBaseline(),parallelVector.normalize()).sum();
}

ParallelPerpendicularBaselineFunctor::OutputType 
ParallelPerpendicularBaselineFunctor
::GetPerpendicularBaseline() const
{
	VectorType perpendicularVector(3);
	perpendicularVector(0) = -sin(m_LookDirection);
	perpendicularVector(1) = 0.0;
	perpendicularVector(2) = cos(m_LookDirection);

	return element_product(this->GetRTNBaseline(),perpendicularVector.normalize()).sum();
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
    case Tangential:
      {
      return this->GetTangentialBaseline();
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
  os << indent << "LookDirection : " << m_LookDirection << std::endl;
}


} // end namespace functor
} // end namespace otb

#endif
