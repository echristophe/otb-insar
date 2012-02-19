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

#include "otbHorizontalVerticalBaselineFunctor.h"


int otbHorizontalVerticalBaselineFunctorTest(int argc, char * argv[])
{

  typedef otb::Functor::HorizontalVerticalBaselineFunctor	BaselineFunctorType;
  typedef BaselineFunctorType::VectorType					VectorType;   

  //Instantiating object
  BaselineFunctorType::Pointer baseline = BaselineFunctorType::New();
  
  VectorType  RTNBaseline(3);
  RTNBaseline[0] = 100.0;
  RTNBaseline[1] = 200.0;
  RTNBaseline[2] = 300.0;

  baseline->SetRTNBaseline(RTNBaseline);
  if(std::abs(baseline->GetHorizontalBaseline() -RTNBaseline[2]) >0.0)
  {
	std::cout << "HorizontalBaseline : " << baseline->GetHorizontalBaseline() << std::endl;
	std::cout << "Expected value : " << RTNBaseline[2] << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetVerticalBaseline() +RTNBaseline[0]) > 0.0 )
  {
	std::cout << "VerticalBaseline : " << baseline->GetVerticalBaseline() << std::endl;
	std::cout << "Expected value : " << -RTNBaseline[0] << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetBaseline(otb::Functor::HorizontalVerticalBaselineFunctor::Horizontal)
	          - baseline->GetHorizontalBaseline()) >0.0)
  {
	std::cout << " Exception when calling GetBaseline() for Horizontal value" << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetBaseline(otb::Functor::HorizontalVerticalBaselineFunctor::Vertical)
	          - baseline->GetVerticalBaseline()) >0.0)
  {
	std::cout << " Exception when calling GetBaseline() for vertical value" << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
