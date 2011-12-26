
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
#include "otbPlatformPositionAdapter.h"

int main(int argc, char* argv[])
{

  if (argc != 3)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile" << std::endl;
    return EXIT_FAILURE;
    }

  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,2> ImageType;

  typedef otb::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(argv[1]);
  slave->SetFileName(argv[2]);

  master->UpdateOutputInformation();
  slave->UpdateOutputInformation();

  typedef otb::PlatformPositionAdapter PlatformType;
  PlatformType::Pointer masterPlatform = PlatformType::New();
  PlatformType::Pointer slavePlatform = PlatformType::New();

  masterPlatform->CreateSensorModel(master->GetOutput()->GetImageKeywordlist());
  slavePlatform->CreateSensorModel(slave->GetOutput()->GetImageKeywordlist());

  std::vector<double> masterPosition;
  std::vector<double> masterSpeed;
  masterPosition.resize(3);
  masterSpeed.resize(3);

  std::vector<double> slavePosition;
  std::vector<double> slaveSpeed;
  slavePosition.resize(3);
  slaveSpeed.resize(3);

  masterPlatform->GetPlatformPosition(0, masterPosition, masterSpeed);
  slavePlatform->GetPlatformPosition(0, slavePosition, slaveSpeed);

  std::cout << std::setprecision(15);

  std::cout << "Master:\n";
  std::cout << masterPosition[0] << "\n";
  std::cout << masterPosition[1] << "\n";
  std::cout << masterPosition[2] << "\n";
  std::cout << masterSpeed[0] << "\n";
  std::cout << masterSpeed[1] << "\n";
  std::cout << masterSpeed[2] << "\n";

  std::cout << "Slave:\n";
  std::cout << slavePosition[0] << "\n";
  std::cout << slavePosition[1] << "\n";
  std::cout << slavePosition[2] << "\n";
  std::cout << slaveSpeed[0] << "\n";
  std::cout << slaveSpeed[1] << "\n";
  std::cout << slaveSpeed[2] << "\n";

  double baselineLength = vcl_sqrt(
      (masterPosition[0] - slavePosition[0]) * (masterPosition[0] - slavePosition[0]) +
      (masterPosition[1] - slavePosition[1]) * (masterPosition[1] - slavePosition[1]) +
      (masterPosition[2] - slavePosition[2]) * (masterPosition[2] - slavePosition[2]));
  std::cout << "Baseline length:\n";
  std::cout << baselineLength << " m \n";

}
