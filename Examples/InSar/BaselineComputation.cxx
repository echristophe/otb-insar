
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

#include <vnl/vnl_vector.h>
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/algo/vnl_lsqr.h>
#include <vnl/vnl_sparse_matrix_linear_system.h>
#include <vnl/vnl_least_squares_function.h>

int main(int argc, char* argv[])
{

  if (argc != 2)
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
  
  typedef otb::PlatformPositionAdapter PlatformType;
  PlatformType::Pointer masterPlatform = PlatformType::New();
  PlatformType::Pointer slavePlatform = PlatformType::New();

  std::vector<double> masterPosition;
  std::vector<double> masterSpeed;
  masterPosition.resize(3);
  masterSpeed.resize(3);

  std::vector<double> slavePosition;
  std::vector<double> slaveSpeed;
  slavePosition.resize(3);
  slaveSpeed.resize(3);
  
  unsigned int numberOfRow  = master->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  unsigned int numberOfCol  = master->GetOutput()->GetLargestPossibleRegion().GetSize()[1];

  std::vector<ImageType::PointType> pointImage;
  pointImage.clear();
  std::vector<double> baselineImage;
  baselineImage.clear();

  for(unsigned int i=0 ; i< numberOfRow ; i+=500)
	{
	for(unsigned int j=0 ; j< numberOfCol ; j+=500)
		{
		ImageType::PointType ImgPoint;

		// References
		ImgPoint[0] = i;
		ImgPoint[1] = j;

		masterPlatform->GetPlatformPosition(ImgPoint[0], masterPosition, masterSpeed);
		slavePlatform->GetPlatformPosition(ImgPoint[0], slavePosition, slaveSpeed);

		std::cout << std::setprecision(15);

		double baselineLength = vcl_sqrt(
					(masterPosition[0] - slavePosition[0]) * (masterPosition[0] - slavePosition[0]) +
					(masterPosition[1] - slavePosition[1]) * (masterPosition[1] - slavePosition[1]) +
					(masterPosition[2] - slavePosition[2]) * (masterPosition[2] - slavePosition[2]));

		//std::cout<<"MasterImgPoint #"<<i<<": ["<<MasterImgPoint[0]<<","<<MasterImgPoint[1]<<"]";
		//std::cout<<" -> estimatedSlaveImgPoint: "<<" ["<<estimatedSlaveImgPoint[0]<<","<<estimatedSlaveImgPoint[1]<<"]";
		//std::cout<<" -> Baseline : " << baselineLength << " m "<<std::endl;
		pointImage.push_back(ImgPoint);
		baselineImage.push_back(baselineLength);
		} 
	}


  std::cout<<"number of points : " << pointImage.size()<<std::endl;

  unsigned int nbPoints = pointImage.size();

  // Convenient typedefs
  typedef vnl_sparse_matrix<double> VnlMatrixType;
  typedef vnl_vector<double>        VnlVectorType;

  /****************************************************************
    // 	Ax = b
    // Balesine(row,col) = a00 + 
    //                     a10*row  + a01*col  +
    //                     a11*row*col +
    //                     a20*row^2 + a02*col^2
   ****************************************************************/

	vnl_sparse_matrix<double> A(nbPoints,6);
	vnl_vector<double> b(nbPoints,0);
	for(unsigned int pId = 0 ; pId < nbPoints ;pId++ )
	{
		A(pId,0) = 1.0;
		A(pId,1) = pointImage.at(pId)[0];
		A(pId,2) = pointImage.at(pId)[1];
		A(pId,3) = pointImage.at(pId)[0] * pointImage.at(pId)[1];
		A(pId,4) = pointImage.at(pId)[0] * pointImage.at(pId)[0];
		A(pId,5) = pointImage.at(pId)[1] * pointImage.at(pId)[1];

		b[pId] = baselineImage[pId];

	}
    // Declare a linear system
    vnl_sparse_matrix_linear_system<double> linearSystem(A, b);

    // A vector where the solution will be stored
    vnl_vector<double> solution(6);

    // Declare the solver
    vnl_lsqr linearSystemSolver(linearSystem);

    // And solve it
    linearSystemSolver.minimize(solution);

	std::cout << "Balesine(row,col) : " << solution[0] 
	          << " + " << solution[1] << " * row"
	          << " + " << solution[2] << " * col"
	          << " + " << solution[3] << " * row*col"
	          << " + " << solution[4] << " * row*row"
	          << " + " << solution[5] << " * col*col"
			  << std::endl;

    double row = 0;
	double col = 0;
    double base =	solution[0] + solution[1]*row + solution[2]*col
					+ solution[3] * row*col
					+ solution[4] * row*row
					+ solution[5] * col*col;    

	std::cout << "(row,col) : " << row << ", " << col << " -> Baseline : " << base << std::endl;

    row = numberOfRow;
	col = numberOfCol;
    base =	solution[0] + solution[1]*row + solution[2]*col
					+ solution[3] * row*col
					+ solution[4] * row*row
					+ solution[5] * col*col;    

	std::cout << "(row,col) : " << row << ", " << col << " -> Baseline : " << base << std::endl;


}
