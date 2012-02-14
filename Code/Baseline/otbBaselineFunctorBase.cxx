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
#ifndef __otbBaselineFunctorBase_cxx
#define __otbBaselineFunctorBase_cxx

#include "otbBaselineFunctorBase.h"

namespace otb
{ 

namespace Functor {


BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetHorizontalBaseline() const
{
	itkExceptionMacro("GetHorizontalBaseline() not implemented");
}

BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetVerticalBaseline() const
{
	itkExceptionMacro("GetVerticalBaseline() not implemented");
}

BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetParallelBaseline() const
{
	itkExceptionMacro("GetParallelBaseline() not implemented");
}

BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetPerpendicularBaseline() const
{
	itkExceptionMacro("GetPerpendicularBaseline() not implemented");
}

BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetLengthBaseline() const
{
	itkExceptionMacro("GetLengthBaseline() not implemented");
}

BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetOrientationBaseline() const
{
	itkExceptionMacro("GetOrientationBaseline() not implemented");
}

BaselineFunctorBase::OutputType 
BaselineFunctorBase
::GetBaseline(BaselineCalculusEnumType map)
{
  switch( map )
    {
    case Horizontal:
      {
      return this->GetHorizontalBaseline();
      break;
      }
    case Vertical:
      {
      return this->GetVerticalBaseline();
      break;
      }
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
    case Length:
      {
      return this->GetLengthBaseline();
      break;
      }
    case Orientation:
      {
      return this->GetOrientationBaseline();
      break;
      }
    }
}

void
BaselineFunctorBase
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << " RTNBaseline vector : " << std::endl;
  os << m_RTNBaseline << std::endl;
}


} // end namespace functor
} // end namespace otb

#endif
