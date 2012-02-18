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
#ifndef __otbPlatformPositionToBaselineCalculator_cxx
#define __otbPlatformPositionToBaselineCalculator_cxx

#include "otbPlatformPositionToBaselineCalculator.h"
#include "otbPlatformPositionAdapter.h"
#include <map>
#include <vnl/vnl_cross.h>
#include "otbHorizontalVerticalBaselineFunctor.h"
#include "otbParallelPerpendicularBaselineFunctor.h"
#include "otbLengthOrientationBaselineFunctor.h"
#include "otbBaselineFunctorBase.h"

namespace otb
{ 
    
/**
 * Constructor
 */
template <class TFunctor>
PlatformPositionToBaselineCalculator<TFunctor>
::PlatformPositionToBaselineCalculator()
{
  m_MasterPlateform = PlatformType::New();
  m_SlavePlateform  = PlatformType::New();
  m_BaselineFunctor = BaselineFunctorType::New();
}


/**
 * Compute Min and Max of m_Image
 */
template <class TFunctor>
typename PlatformPositionToBaselineCalculator<TFunctor>::BaselineFunctorOutputType
PlatformPositionToBaselineCalculator<TFunctor>
::Evaluate(double line,BaselineCalculusEnumType map) 
{
	std::vector<double> masterPosition(3);
	std::vector<double> slavePosition(3);
	std::vector<double> masterSpeed(3);
	std::vector<double> slaveSpeed(3);

	m_MasterPlateform->GetPlatformPosition(line, masterPosition, masterSpeed);
	m_SlavePlateform->GetPlatformPosition(line, slavePosition, slaveSpeed);

	vnl_vector<double> baselineVector(3);

	baselineVector = this->BaselineInRTNSystem(masterPosition, slavePosition, masterSpeed);
	m_BaselineFunctor->SetRTNBaseline(baselineVector);
	return this->m_BaselineFunctor->GetBaseline(map);
}




/**
 * Evaluate Baseline in RTN (Radial Tangential Normal) System coordinate 
 */
template <class TFunctor>
vnl_vector<double>
PlatformPositionToBaselineCalculator<TFunctor>
::BaselineInRTNSystem(
				std::vector<double> & masterPosition,
				std::vector<double> & slavePosition,
				std::vector<double> & masterSpeed)
{
	vnl_vector<double> baselineRTN(3);

	/** Define the Radial vector */
	vnl_vector<double> radialVector(3);
	radialVector(0) = masterPosition[0];
	radialVector(1) = masterPosition[1];
	radialVector(2) = masterPosition[2];
	radialVector.normalize();

    /** Define the Tagential vector */
	vnl_vector<double> tangentialVector(3);
	tangentialVector(0) = masterSpeed[0];
	tangentialVector(1) = masterSpeed[1];
	tangentialVector(2) = masterSpeed[2];
	tangentialVector.normalize();

    /** Define the Normal vector */
	vnl_vector<double> normalVector(3);
	normalVector = vnl_cross_3d(radialVector,tangentialVector);

	vnl_vector<double> baselineXYZ(3);

	baselineXYZ(0) = masterPosition[0] - slavePosition[0];
	baselineXYZ(1) = masterPosition[1] - slavePosition[1];
	baselineXYZ(2) = masterPosition[2] - slavePosition[2];

	//baselineRTN
	baselineRTN[0] = dot_product(baselineXYZ,radialVector); 
	baselineRTN[1] = dot_product(baselineXYZ,tangentialVector); 
	baselineRTN[2] = dot_product(baselineXYZ,normalVector); 

	return baselineRTN;
}


template <class TFunctor>
void
PlatformPositionToBaselineCalculator<TFunctor>
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "MasterPlateform: " << std::endl;
  m_MasterPlateform->Print(os, indent.GetNextIndent());
  os << indent << "m_SlavePlateform: " << std::endl;
  m_SlavePlateform->Print(os, indent.GetNextIndent());

}

} // end namespace otb

#endif
