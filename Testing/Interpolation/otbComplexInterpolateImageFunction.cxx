/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See OTBCopyright.txt for details.


  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkMacro.h"

#include "otbComplexInterpolateImageFunction.h"
#include "otbImage.h"
#include "otbImageFileReader.h"
#include <fstream>

namespace Function
{
template<class TInput = double, class TOutput = double>
class ComplexFunction
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

int otbComplexInterpolateImageFunction(int argc, char * argv[])
{
  const char * infname = argv[1];
  const char * outfname = argv[2];

  typedef std::complex<double>										InputPixelType;
  typedef otb::Image<InputPixelType, 2>								ImageType;
  typedef Function::ComplexFunction<InputPixelType, InputPixelType>	FunctionType;
  typedef itk::ConstantBoundaryCondition<ImageType>					BoundaryConditionType;
  typedef double													CoordRepType;

  typedef otb::ComplexInterpolateImageFunction<ImageType, 
						FunctionType, BoundaryConditionType,
						CoordRepType>								InterpolatorType;

  typedef InterpolatorType::ContinuousIndexType                     ContinuousIndexType;
  typedef otb::ImageFileReader<ImageType>                           ReaderType;

  unsigned int radius = std::atoi(argv[3]);
  double zeroFrequency = std::atof(argv[4]);

  int i = 5;

  std::vector<ContinuousIndexType> indicesList;

  while (i < argc && (i + 1) < argc)
    {
    ContinuousIndexType idx;

    idx[0] = atof(argv[i]);
    idx[1] = atof(argv[i + 1]);

    indicesList.push_back(idx);

    i += 2;
    }

  // Instantiating object
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(infname);
  reader->Update();

  interpolator->SetInputImage(reader->GetOutput());
  interpolator->SetRadius(radius);
  interpolator->SetZeroFrequency(zeroFrequency);

  std::ofstream file;
  file.open(outfname);

  for (std::vector<ContinuousIndexType>::iterator it = indicesList.begin(); it != indicesList.end(); ++it)
    {
    file << (*it) << " -> " << interpolator->EvaluateAtContinuousIndex((*it)) << std::endl;
    }

  file.close();

  return EXIT_SUCCESS;
}
