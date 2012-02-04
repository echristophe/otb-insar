
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
#include "otbGenericRSTransform.h"


int main(int argc, char* argv[])
{

  if (argc != 3)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile demDir" << std::endl;
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
  
//  std::string demDir    = argv[3];

  // Build wgs ref
  OGRSpatialReference oSRS;
  oSRS.SetWellKnownGeogCS("WGS84");
  char * wgsRef = NULL;
  oSRS.exportToWkt(&wgsRef);

  typedef otb::GenericRSTransform<>         TransformType;

  // Instantiate MasterImage->WGS transform
  TransformType::Pointer masterImg2wgs = TransformType::New();
  masterImg2wgs->SetInputProjectionRef(master->GetOutput()->GetProjectionRef());
  masterImg2wgs->SetInputKeywordList(master->GetOutput()->GetImageKeywordlist());
  masterImg2wgs->SetOutputProjectionRef(wgsRef);
//  masterImg2wgs->SetDEMDirectory(demDir);
  masterImg2wgs->InstanciateTransform();

  // Instantiate WGS->Image transform
  TransformType::Pointer wgs2SlaveImg = TransformType::New();
  wgs2SlaveImg->SetInputProjectionRef(wgsRef);
  wgs2SlaveImg->SetOutputProjectionRef(slave->GetOutput()->GetProjectionRef());
  wgs2SlaveImg->SetOutputKeywordList(slave->GetOutput()->GetImageKeywordlist());
//  wgs2SlaveImg->SetDEMDirectory(demDir);
  wgs2SlaveImg->InstanciateTransform();



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
  
  unsigned int numberOfLine  = master->GetOutput()->GetLargestPossibleRegion().GetSize()[0];


  for(unsigned int i=0 ; i< numberOfLine ; i+=100)
  {
	  ImageType::PointType MasterImgPoint, estimatedMasterGeoPoint, estimatedSlaveImgPoint;

    // References
    MasterImgPoint[0] = i;
    MasterImgPoint[1] = 0;

	// Estimations Slave image point
    estimatedMasterGeoPoint = masterImg2wgs->TransformPoint(MasterImgPoint);
    estimatedSlaveImgPoint = wgs2SlaveImg->TransformPoint(estimatedMasterGeoPoint);


	if(estimatedSlaveImgPoint[0] <0)
	{
		std::cout<<"MasterImgPoint #"<<i<<": ["<<MasterImgPoint[0]<<","<<MasterImgPoint[1]<<"]";
		std::cout<<" -> estimatedSlaveImgPoint: "<<" ["<<estimatedSlaveImgPoint[0]<<","<<estimatedSlaveImgPoint[1]<<"]";
		std::cout<<" -> No Baseline "<<std::endl;
	}
	else
	{
		masterPlatform->GetPlatformPosition(MasterImgPoint[0], masterPosition, masterSpeed);
		slavePlatform->GetPlatformPosition(estimatedSlaveImgPoint[0], slavePosition, slaveSpeed);

		std::cout << std::setprecision(15);

	double baselineLength = vcl_sqrt(
			(masterPosition[0] - slavePosition[0]) * (masterPosition[0] - slavePosition[0]) +
			(masterPosition[1] - slavePosition[1]) * (masterPosition[1] - slavePosition[1]) +
			(masterPosition[2] - slavePosition[2]) * (masterPosition[2] - slavePosition[2]));

	std::cout<<"MasterImgPoint #"<<i<<": ["<<MasterImgPoint[0]<<","<<MasterImgPoint[1]<<"]";
	std::cout<<" -> estimatedSlaveImgPoint: "<<" ["<<estimatedSlaveImgPoint[0]<<","<<estimatedSlaveImgPoint[1]<<"]";
	std::cout<<" -> Baseline : " << baselineLength << " m "<<std::endl;
	}
  }

}
