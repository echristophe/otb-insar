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

#include "otbImage.h"
#include "otbBaselineCalculator.h"
#include "otbPlatformPositionToBaselineCalculator.h"
#include "otbBaselineFunctorBase.h"

#include <vnl/vnl_sparse_matrix.h>

namespace otb
{ 
    
/**
 * Constructor
 */
template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
BaselineCalculator<TMasterInputImage,TSlaveInputImage, TFunctor>
::BaselineCalculator()
{
  m_MasterImage = TMasterInputImage::New();
  m_SlaveImage = TSlaveInputImage::New();
  m_PlateformPositionToBaselineCalculator = BaselineType::New();
  m_BaselineCoefficient.clear();
}


/**
 * Compute Min and Max of m_Image
 */
template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
void
BaselineCalculator<TMasterInputImage,TSlaveInputImage, TFunctor>
::Compute(BaselineCalculusEnumType map) 
{
  std::vector<typename MasterImageType::PointType> pointImage;
  pointImage.clear();
  std::vector<double> baselineImage;
  baselineImage.clear();

  this->ExtractBaseline(map, pointImage, baselineImage);

  m_BaselineCoefficient = this->BaselineLinearSolve(pointImage,baselineImage);
}


template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
void
BaselineCalculator<TMasterInputImage,TSlaveInputImage, TFunctor>
::ExtractBaseline(  BaselineCalculusEnumType map,
					std::vector<typename MasterImageType::PointType> & pointImage,
					std::vector<double> & baselineImage) 
{
  unsigned int numberOfRow  = m_MasterImage->GetLargestPossibleRegion().GetSize()[0];
  unsigned int numberOfCol  = m_MasterImage->GetLargestPossibleRegion().GetSize()[1];

  typename BaselineType::Pointer baselineCalculator = BaselineType::New();
  baselineCalculator->SetMasterPlateform(m_MasterImage->GetImageKeywordlist());
  baselineCalculator->SetSlavePlateform(m_SlaveImage->GetImageKeywordlist());

  std::vector<double> lengthBaselineImage;
  lengthBaselineImage.clear();

  for(unsigned int i=0 ; i< numberOfRow ; i+=500)
	{
	for(unsigned int j=0 ; j< numberOfCol ; j+=500)
		{
		typename MasterImageType::PointType ImgPoint;

		// References
		ImgPoint[0] = i;
		ImgPoint[1] = j;

		double baselineValue = baselineCalculator->Evaluate(ImgPoint[0],map);
		pointImage.push_back(ImgPoint);
		baselineImage.push_back(baselineValue);
		}
	}

  std::cout<<"number of points : " << pointImage.size()<<std::endl;

}


template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
vnl_vector<double>
BaselineCalculator<TMasterInputImage,TSlaveInputImage, TFunctor>
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

	vnl_matrix<double> A(nbPoints,6);
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

    // A vector where the solution will be stored
    vnl_vector<double> solution(6);
	solution = vnl_matrix_inverse<double>(A.transpose()*A)*(A.transpose()*b);

	return solution;
}

template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
double
BaselineCalculator<TMasterInputImage,TSlaveInputImage, TFunctor>
::EvaluateBaseline(double row,double col)
{
	double result = 0;
	result =	  m_BaselineCoefficient.get(0) 
				+ m_BaselineCoefficient.get(1) * row 
				+ m_BaselineCoefficient.get(2) * col
				+ m_BaselineCoefficient.get(3) * row*col
				+ m_BaselineCoefficient.get(4) * row*row
				+ m_BaselineCoefficient.get(5) * col*col;   

	return result;
}


template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
void
BaselineCalculator<TMasterInputImage,TSlaveInputImage, TFunctor>
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "MasterImage: " << std::endl;
  m_MasterImage->Print(os, indent.GetNextIndent());
  os << indent << "SlaveImage: " << std::endl;
  m_SlaveImage->Print(os, indent.GetNextIndent());
  os << indent << "SlaveImage: " << m_BaselineCoefficient << std::endl;

}

} // end namespace otb

#endif
