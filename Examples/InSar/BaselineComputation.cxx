
/*=========================================================================

   Copyright 2011 Emmanuel Christophe
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

// Command line:
//
// ./bin/BaselineComputation ~/project/Images/TSX1_SAR__SSC______HS_S_SRA_20090212T204239_20090212T204240/TSX1_SAR__SSC______HS_S_SRA_20090212T204239_20090212T204240.xml ~/project/Images/TSX1_SAR__SSC______HS_S_SRA_20090223T204240_20090223T204241/TSX1_SAR__SSC______HS_S_SRA_20090223T204240_20090223T204241.xml

#include <iomanip>

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbBaselineCalculator.h"
#include "otbLengthOrientationBaselineFunctor.h"
#include "otbPlatformPositionToBaselineCalculator.h"

int main(int argc, char* argv[])
{

  if (argc != 3)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 2 ;


  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,Dimension> ImageType;

  typedef otb::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(argv[1]);
  slave->SetFileName(argv[2]);

  master->UpdateOutputInformation();
  slave->UpdateOutputInformation();
  typedef otb::Functor::LengthOrientationBaselineFunctor	BaselineFunctorType;
  typedef otb::BaselineCalculator<BaselineFunctorType>    BaselineCalculatorType;
  typedef BaselineCalculatorType::PlateformPositionToBaselineCalculatorType PlateformPositionToBaselineCalculatorType;
  BaselineCalculatorType::Pointer baselineCalculator = BaselineCalculatorType::New();

  BaselineCalculatorType::PlateformPositionToBaselinePointer plateformPositionToBaseline =  PlateformPositionToBaselineCalculatorType::New();
  plateformPositionToBaseline->SetMasterPlateform(master->GetOutput()->GetImageKeywordlist());
  plateformPositionToBaseline->SetSlavePlateform(slave->GetOutput()->GetImageKeywordlist());

  baselineCalculator->SetPlateformPositionToBaselineCalculator(plateformPositionToBaseline);

  baselineCalculator->Compute(otb::Functor::LengthOrientationBaselineFunctor::Length);

  double row = 0;
  double col = 0;
  std::cout << "(row,col) : " << row << ", " << col << " -> Baseline : ";
  std::cout <<  baselineCalculator->EvaluateBaseline(row,col)<< std::endl;

  row = 1000;
  col = 1000;
  std::cout << "(row,col) : " << row << ", " << col << " -> Baseline : ";
  std::cout <<  baselineCalculator->EvaluateBaseline(row,col)<< std::endl;
}
