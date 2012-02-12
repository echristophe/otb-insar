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
#ifndef __otbLengthOrientationBaselineFunctor_cxx
#define __otbLengthOrientationBaselineFunctor_cxx

#include "otbLengthOrientationBaselineFunctor.h"

namespace otb
{ 

namespace Functor {


  
LengthOrientationBaselineFunctor::MapType 
LengthOrientationBaselineFunctor
::operator ()(const VectorPositionType & baselineRTN) const
{
    MapType out;
	out.clear();

	double baselineLength = baselineRTN.two_norm(); 
  	out.insert(std::pair<std::string,double>("Length",baselineLength) );	
	vnl_vector<double> normalComponent(3);
	normalComponent.fill(0.0);
	normalComponent(2) = 1.0;
	double angle = acos(dot_product(baselineRTN,normalComponent) / baselineLength) * CONST_180_PI; 
  	out.insert(std::pair<std::string,double>("Angle",angle) );	
	return out;
}

LengthOrientationBaselineFunctor::MapType 
LengthOrientationBaselineFunctor
::operator ()(const vnl_vector<double> & val1,const vnl_vector<double> & val2) const
{
	itkExceptionMacro("Only one argument to operator()");
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
