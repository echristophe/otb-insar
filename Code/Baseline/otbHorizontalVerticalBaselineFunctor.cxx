/*=========================================================================

   Copyright 2012 Patrick IMBO
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
#ifndef __otbHorizontalVerticalBaselineFunctor_cxx
#define __otbHorizontalVerticalBaselineFunctor_cxx

#include "otbHorizontalVerticalBaselineFunctor.h"
#include "otbBaselineFunctorBase.h"

namespace otb
{ 

namespace Functor {

HorizontalVerticalBaselineFunctor::OutputType
HorizontalVerticalBaselineFunctor
::GetHorizontalBaseline() const
{
	OutputType horizontalBaseline;
	horizontalBaseline = this->GetRTNBaseline().get(2);
	return horizontalBaseline;
}


HorizontalVerticalBaselineFunctor::OutputType
HorizontalVerticalBaselineFunctor
::GetVerticalBaseline() const
{
	OutputType verticalBaseline;
	verticalBaseline = this->GetRTNBaseline().get(1);
	return verticalBaseline;
}


void
HorizontalVerticalBaselineFunctor
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}


} // end namespace functor
} // end namespace otb

#endif
