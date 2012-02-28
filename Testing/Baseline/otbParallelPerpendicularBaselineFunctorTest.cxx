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

#include "otbParallelPerpendicularBaselineFunctor.h"
#include "otbMath.h"
#include "vnl/vnl_cross.h"

int otbParallelPerpendicularBaselineFunctorTest(int argc, char * argv[])
{

  typedef otb::Functor::ParallelPerpendicularBaselineFunctor	BaselineFunctorType;
  typedef BaselineFunctorType::VectorType					VectorType;   

  //Instantiating object
  BaselineFunctorType::Pointer baseline = BaselineFunctorType::New();
  
  VectorType  RTNBaseline(3);
  RTNBaseline.fill(0.0);

  const double angle = 0.0; // in degree
  const double parallelBaseline = 100.0;
  const double perpendicularBaseline = 52.0;
  double degToRad = angle *otb::CONST_PI_180;

  VectorType  refPoint(3);
  degToRad = angle *otb::CONST_PI_180;
  refPoint[0] = parallelBaseline * cos(degToRad);
  refPoint[1] = 0.0;
  refPoint[2] = parallelBaseline * sin(degToRad);

  std::cout << "refPoint : " << refPoint << std::endl;

  VectorType  tangentialVect(3);
  tangentialVect.fill(0.0);
  tangentialVect(1) = 1.0;

  std::cout << "tangentialVect : " << tangentialVect << std::endl;

  VectorType  perpendicularVect(3);
  perpendicularVect = perpendicularBaseline * vnl_cross_3d(refPoint,tangentialVect).normalize();

  std::cout << "perpendicularVect : " << perpendicularVect << std::endl;

  RTNBaseline = refPoint + perpendicularVect;

  std::cout << "RTNBaseline : " << RTNBaseline << std::endl;

  baseline->SetRTNBaseline(RTNBaseline);
  baseline->SetRefPoint(refPoint);

  if(std::abs(baseline->GetParallelBaseline() - parallelBaseline) > 0.000001 )
  {
	std::cout << "parallelBaseline : " << baseline->GetParallelBaseline() << std::endl;
	std::cout << "Expected value : " << parallelBaseline << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetPerpendicularBaseline() -perpendicularBaseline) > 0.0)
  {
	std::cout << "perpendicularBaselineBaseline : " << baseline->GetPerpendicularBaseline() << std::endl;
	std::cout << "Expected value : " << perpendicularBaseline << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetBaseline(otb::Functor::ParallelPerpendicularBaselineFunctor::Parallel)
				- baseline->GetParallelBaseline()) >0.0)
  {
	std::cout << " Exception when calling GetBaseline() for Parallel enum" << std::endl;
    return EXIT_FAILURE;
  }

  if(std::abs(baseline->GetBaseline(otb::Functor::ParallelPerpendicularBaselineFunctor::Perpendicular)
				- baseline->GetPerpendicularBaseline()) >0.0)
  {
	std::cout << " Exception when calling GetBaseline() for Perpendicular enum" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
