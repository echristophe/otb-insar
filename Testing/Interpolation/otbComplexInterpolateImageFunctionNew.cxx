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
#include "itkMacro.h"

#include "otbComplexInterpolateImageFunction.h"
#include "otbImage.h"
#include "itkConstantBoundaryCondition.h"

namespace Function
{
template<class TInput = double, class TOutput = double>
class SameFunction
{
public:
  void SetRadius(unsigned int rad)
  {
    m_Radius = rad;
  }
  unsigned int GetRadius() const
  {
    return m_Radius;
  }
  inline TOutput operator ()(const TInput& A) const
  {
    return static_cast<TOutput>(A);
  }
  unsigned int m_Radius;
};

}

int otbComplexInterpolateImageFunctionNew(int argc, char * argv[])
{

  typedef std::complex<double> InputPixelType;
  const int Dimension = 2;
  typedef otb::Image<InputPixelType, Dimension>                  ImageType;
  typedef Function::SameFunction<InputPixelType, InputPixelType> FunctionType;
  typedef itk::ConstantBoundaryCondition<ImageType>              BoundaryConditionType;
  typedef double                                                 CoordRepType;

  typedef otb::ComplexInterpolateImageFunction<ImageType, FunctionType, BoundaryConditionType,
      CoordRepType> GenericFunctionType;

  // Instantiating object
  GenericFunctionType::Pointer generic = GenericFunctionType::New();

  std::cout << generic << std::endl;

  return EXIT_SUCCESS;
}
