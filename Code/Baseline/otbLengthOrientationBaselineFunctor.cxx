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
#ifndef __otbLengthOrientationBaselineFunctor_cxx
#define __otbLengthOrientationBaselineFunctor_cxx

#include "otbLengthOrientationBaselineFunctor.h"
#include "otbMath.h"

namespace otb
{ 

namespace Functor {



LengthOrientationBaselineFunctor::OutputType
LengthOrientationBaselineFunctor
::GetLengthBaseline() const
{
	double baselineLength = this->GetRTNBaseline().two_norm(); 
	return baselineLength;
}

LengthOrientationBaselineFunctor::OutputType
LengthOrientationBaselineFunctor
::GetOrientationBaseline() const
{
	vnl_vector<double> normalComponent(3);
	normalComponent.fill(0.0);
	normalComponent(2) = 1.0;
	double baselineLength = this->GetLengthBaseline();
	double angle = acos(dot_product(this->GetRTNBaseline(),normalComponent) / baselineLength) * CONST_180_PI; 
	return angle;
}

LengthOrientationBaselineFunctor::OutputType
LengthOrientationBaselineFunctor
::GetBaseline(BaselineCalculusEnumType map)
{
  switch( map )
    {
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
	return 0;
}


void
LengthOrientationBaselineFunctor
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}


} // end namespace functor
} // end namespace otb

#endif
