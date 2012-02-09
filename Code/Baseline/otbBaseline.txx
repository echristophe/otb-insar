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

namespace otb
{ 
    
/**
 * Constructor
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::Baseline()
{
  m_MasterImage = TMasterInputImage::New();
  m_SlaveImage = TSlaveInputImage::New();
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
	m_Baseline = m_Functor(masterPosition,slavePosition);
}



/**
 * Compute Master plateform position
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
std::vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::GetMasterPlateformPosition(double line)
{
	typedef otb::PlatformPositionAdapter PlatformType;
	PlatformType::Pointer platform = PlatformType::New();
	
	platform->CreateSensorModel(m_MasterImage->GetImageKeywordlist());

    std::vector<double> position;
	std::vector<double> speed;
	position.resize(3);
	speed.resize(3);

	platform->GetPlatformPosition(line, position, speed);

	return position;
}

/**
 * Compute Slave plateform position
 */
template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
std::vector<double>
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
::GetSlavePlateformPosition(double line)
{
	typedef otb::PlatformPositionAdapter PlatformType;
	PlatformType::Pointer platform = PlatformType::New();
	
	platform->CreateSensorModel(m_SlaveImage->GetImageKeywordlist());

	std::vector<double> position;
	std::vector<double> speed;
	position.resize(3);
	speed.resize(3);

	platform->GetPlatformPosition(line, position, speed);

	return position;
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


template <class TMasterInputImage,class TSlaveInputImage,class TBaselineFunctor>
void
Baseline<TMasterInputImage,TSlaveInputImage,TBaselineFunctor>
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
