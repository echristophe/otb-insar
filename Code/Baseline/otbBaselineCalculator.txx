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
#ifndef __otbBaselineCalculator_txx
#define __otbBaselineCalculator_txx

#include "otbBaselineCalculator.h"
#include "otbBaseline.h"
#include "otbImage.h"

#include <vnl/vnl_sparse_matrix.h>
#include <vnl/algo/vnl_lsqr.h>
#include <vnl/vnl_sparse_matrix_linear_system.h>
#include <vnl/vnl_least_squares_function.h>

namespace otb
{ 
    
/**
 * Constructor
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
BaselineCalculator<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::BaselineCalculator()
{
  m_MasterImage = TMasterInputImage::New();
  m_SlaveImage = TSlaveInputImage::New();
  m_BaselineCoefficient.clear();
}


/**
 * Compute Min and Max of m_Image
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void
BaselineCalculator<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::Compute(void) 
{
  std::vector<MasterImageType::PointType> pointImage;
  pointImage.clear();
  std::map<std::string,std::vector<double> > baselineImage;
  baselineImage.clear();

  this->ExtractBaseline( pointImage, baselineImage);

  std::map<std::string,std::vector<double> >::const_iterator it;
  it = baselineImage.find("Length");
  std::vector<double> value = it->second; 

  vnl_vector<double> coefLength;
  coefLength = this->BaselineLinearSolve(pointImage,value);
  m_BaselineCoefficient["Length"]= coefLength;
}


template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void
BaselineCalculator<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::ExtractBaseline(std::vector<typename MasterImageType::PointType> & pointImage,
				  std::map<std::string,std::vector<double> > & baselineImage) 
{
  unsigned int numberOfRow  = m_MasterImage->GetLargestPossibleRegion().GetSize()[0];
  unsigned int numberOfCol  = m_MasterImage->GetLargestPossibleRegion().GetSize()[1];

  BaselineType::Pointer baselineCalculator = BaselineType::New();
  baselineCalculator->SetMasterPlateform(m_MasterImage->GetImageKeywordlist());
  baselineCalculator->SetSlavePlateform(m_SlaveImage->GetImageKeywordlist());

  std::vector<double> lengthBaselineImage;
  lengthBaselineImage.clear();

  for(unsigned int i=0 ; i< numberOfRow ; i+=500)
	{
	for(unsigned int j=0 ; j< numberOfCol ; j+=500)
		{
		MasterImageType::PointType ImgPoint;

		// References
		ImgPoint[0] = i;
		ImgPoint[1] = j;

		baselineCalculator->Compute(ImgPoint[0]);
		double baselineLength = baselineCalculator->GetBaselineValue("Length"); 
		lengthBaselineImage.push_back(baselineLength);
		pointImage.push_back(ImgPoint);
		//baselineImage.push_back(baselineLength);
		}
	}
  baselineImage["Length"] = lengthBaselineImage;

  std::cout<<"number of points : " << pointImage.size()<<std::endl;

}


template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
vnl_vector<double>
BaselineCalculator<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::BaselineLinearSolve(std::vector<typename MasterImageType::PointType> & pointImage,
		   std::vector<double> & baselineImage) 
{
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

	return solution;
}

template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
double
BaselineCalculator<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::EvaluateBaseline(double row,double col)
{
    vnl_vector<double> solution(6);
	CoefMapType::const_iterator it;

	double result = 0;
	for(it = m_BaselineCoefficient.begin(); it != m_BaselineCoefficient.end(); ++it)
	{
	std::cout << "Key: " << (*it).first; //<< " Value: " << (*itr).second;
	//it = m_BaselineCoefficient.find("Length");
	solution = it->second; 

	result =	solution[0] + solution[1]*row + solution[2]*col
					+ solution[3] * row*col
					+ solution[4] * row*row
					+ solution[5] * col*col;   
	}
	return result;
}


template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void
BaselineCalculator<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "MasterImage: " << std::endl;
  m_MasterImage->Print(os, indent.GetNextIndent());
  os << indent << "SlaveImage: " << std::endl;
  m_SlaveImage->Print(os, indent.GetNextIndent());

}

} // end namespace otb

#endif
