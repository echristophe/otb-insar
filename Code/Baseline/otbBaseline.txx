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
#ifndef __otbBaseline_txx
#define __otbBaseline_txx

#include "otbBaseline.h"
#include "otbPlatformPositionAdapter.h"
#include <map>
#include <vnl/vnl_cross.h>

namespace otb
{ 
    
/**
 * Constructor
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::Baseline()
{
  m_MasterPlateform = PlatformType::New();
  m_SlavePlateform  = PlatformType::New();
  m_Baseline.clear();
}


/**
 * Compute Min and Max of m_Image
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::Compute(double line) 
{
	std::vector<double> masterPosition;
	std::vector<double> slavePosition;
	masterPosition.resize(3);
	slavePosition.resize(3);

	this->EvaluateMasterAndSlavePosition(line,line, masterPosition, slavePosition);

	std::vector<double> masterSpeed;
	std::vector<double> slaveSpeed;
	masterSpeed.resize(3);
	slaveSpeed.resize(3);

	this->EvaluateMasterAndSlaveSpeed(line,line, masterSpeed, slaveSpeed);

	vnl_vector<double> baselineVector(3);

	baselineVector = this->BaselineInRTNSystem(masterPosition, slavePosition, masterSpeed);
	m_Baseline = m_Functor(baselineVector);
}



/**
 * Compute Master plateform position
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
std::vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::GetMasterPlateformPosition(double line)
{
    std::vector<double> position;
	std::vector<double> speed;
	position.resize(3);
	speed.resize(3);

	m_MasterPlateform->GetPlatformPosition(line, position, speed);

	return position;
}

/**
 * Compute Master plateform position
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
std::vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::GetMasterPlateformSpeed(double line)
{
    std::vector<double> position;
	std::vector<double> speed;
	position.resize(3);
	speed.resize(3);

	m_MasterPlateform->GetPlatformPosition(line, position, speed);

	return speed;
}

/**
 * Compute Slave plateform position
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
std::vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::GetSlavePlateformPosition(double line)
{
	std::vector<double> position;
	std::vector<double> speed;
	position.resize(3);
	speed.resize(3);

	m_SlavePlateform->GetPlatformPosition(line, position, speed);

	return position;
}

/**
 * Compute Slave plateform speed
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
std::vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::GetSlavePlateformSpeed(double line)
{
	std::vector<double> position;
	std::vector<double> speed;
	position.resize(3);
	speed.resize(3);

	m_SlavePlateform->GetPlatformPosition(line, position, speed);

	return speed;
}

/**
 * Evaluate Master and Slave plateform position
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void 
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::EvaluateMasterAndSlavePosition(
				double masterLine, double slaveLine,
				std::vector<double> & masterPosition,
				std::vector<double> & slavePosition)
{
    masterPosition = this->GetMasterPlateformPosition(masterLine);
    slavePosition  = this->GetSlavePlateformPosition(slaveLine);
}

/**
 * Evaluate Master and Slave plateform speed
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void 
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::EvaluateMasterAndSlaveSpeed(
				double masterLine, double slaveLine,
				std::vector<double> & masterSpeed,
				std::vector<double> & slaveSpeed)
{
    masterSpeed = this->GetMasterPlateformSpeed(masterLine);
    slaveSpeed  = this->GetSlavePlateformSpeed(slaveLine);
}


/**
 * Evaluate Baseline in RTN (Radial Tangential Normal) System coordinate 
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
vnl_vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
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



template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
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
