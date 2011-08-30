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
#include "otbImageFileReader.h"
#include <fstream>

#include "otbWindowedSincInterpolateImageBlackmanFunction.h"
#include "otbWindowedSincInterpolateImageHammingFunction.h"

template<class TFunction>
int otbComplexInterpolateImageFunction_generic(int argc, char * argv[])
{
  const char * infname = argv[1];
  const char * outfname = argv[2];

  typedef std::complex<double>                    InputPixelType;
  typedef otb::Image<InputPixelType, 2>                ImageType;
  typedef TFunction                            FunctionType;
  typedef itk::ConstantBoundaryCondition<ImageType>          BoundaryConditionType;
  typedef double                          CoordRepType;

  typedef otb::ComplexInterpolateImageFunction<ImageType, 
            FunctionType, BoundaryConditionType,
            CoordRepType>                InterpolatorType;

  typedef InterpolatorType::ContinuousIndexType                     ContinuousIndexType;
  typedef otb::ImageFileReader<ImageType>                           ReaderType;

  unsigned int radius = std::atoi(argv[4]);
  double zeroFrequency = std::atof(argv[5]);

  int i = 6;

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
  interpolator->SetNormalizeZeroFrequency(zeroFrequency);

  std::ofstream file;
  file.open(outfname);

  for (std::vector<ContinuousIndexType>::iterator it = indicesList.begin(); it != indicesList.end(); ++it)
    {
    file << (*it) << " -> " << interpolator->EvaluateAtContinuousIndex((*it)) << std::endl;
    }

  file.close();

  return EXIT_SUCCESS;
}


int otbComplexInterpolateImageFunction(int argc, char * argv[])
{
  int FunctionType = atoi(argv[3]);
  switch (FunctionType)
    {
    case 0:
      return otbComplexInterpolateImageFunction_generic<otb::Function::BlackmanWindowFunction<double> > (argc, argv);
      break;
    case 1:
      return otbComplexInterpolateImageFunction_generic<otb::Function::HammingWindowFunction<double> > (argc, argv);
      break;
    default:
      std::cerr << "No more function available\n";
      return EXIT_FAILURE;
    }
}