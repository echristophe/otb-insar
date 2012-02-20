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
#include "itkMacro.h"

#include "otbLengthOrientationBaselineFunctor.h"
#include "otbMath.h"


int otbLengthOrientationBaselineFunctorTest(int argc, char * argv[])
{

  typedef otb::Functor::LengthOrientationBaselineFunctor	BaselineFunctorType;
  typedef BaselineFunctorType::VectorType					VectorType;   

  //Instantiating object
  BaselineFunctorType::Pointer baseline = BaselineFunctorType::New();
  
  VectorType  RTNBaseline(3);
  double angle = 35.0; // in degree
  double length = 100.0;
  double degToRad = angle *otb::CONST_PI_180;
  RTNBaseline[0] = -length * sin(degToRad);
  RTNBaseline[1] = 0.0;
  RTNBaseline[2] = length * cos(degToRad);

  baseline->SetRTNBaseline(RTNBaseline);


  double resultLength = std::sqrt(	  RTNBaseline[0]*RTNBaseline[0]
									+ RTNBaseline[1]*RTNBaseline[1]
									+ RTNBaseline[2]*RTNBaseline[2]);

  if(std::abs(baseline->GetLengthBaseline() - length) > 0.000001 )
  {
	std::cout << "LengthBaseline : " << baseline->GetLengthBaseline() << std::endl;
	std::cout << "Expected value : " << length << std::endl;
    return EXIT_FAILURE;
  }

  double resultOrientation = std::atan2(-RTNBaseline[0],-RTNBaseline[2])* otb::CONST_180_PI;

  if(std::abs(baseline->GetOrientationBaseline() -angle) > 0.0)
  {
	std::cout << "OrientationBaseline : " << baseline->GetOrientationBaseline() << std::endl;
	std::cout << "Expected value : " << angle << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetBaseline(otb::Functor::LengthOrientationBaselineFunctor::Orientation)
				- baseline->GetOrientationBaseline()) >0.0)
  {
	std::cout << " Exception when calling GetBaseline() for Orientation enum" << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetBaseline(otb::Functor::LengthOrientationBaselineFunctor::Length)
				- baseline->GetLengthBaseline()) >0.0)
  {
	std::cout << " Exception when calling GetBaseline() for Length enum" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
