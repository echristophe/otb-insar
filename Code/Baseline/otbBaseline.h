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
#ifndef __otbBaseline_h
#define __otbBaseline_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "otbMacro.h"
#include <vnl/vnl_vector.h>
#include "otbPlatformPositionAdapter.h"

namespace otb
{

namespace Functor {

class LengthOrientationBaseline
{
public:
  typedef std::vector<double> PositionType;
  typedef std::map<std::string,double> MapType;
  LengthOrientationBaseline() {}
  virtual ~LengthOrientationBaseline() {}

  inline MapType operator ()(const vnl_vector<double> & baselineRTN)
  {
    MapType out;
	out.clear();

	double baselineLength = baselineRTN.two_norm(); 
  	out.insert(std::pair<std::string,double>("Length",baselineLength) );	
	vnl_vector<double> normalComponent(3);
	normalComponent.fill(0.0);
	normalComponent(2) = 1.0;
	double angle = acos(dot_product(baselineRTN,normalComponent) / baselineLength) * CONST_180_PI; 
  	out.insert(std::pair<std::string,double>("Angle",angle) );	
	return out;
  }
};

class ParallelPerpendicularBaseline
{
public:
  typedef std::vector<double> PositionType;
  typedef std::map<std::string,double> MapType;
  ParallelPerpendicularBaseline() {}
  virtual ~ParallelPerpendicularBaseline() {}


  inline MapType operator ()(const vnl_vector<double> & baselineRTN)
  {
    MapType out;
	out.clear();
	/** TODO*/

	return out;
  }
};

class HorizontalVerticalBaseline
{
public:
  typedef std::vector<double> PositionType;
  typedef std::map<std::string,double> MapType;
  HorizontalVerticalBaseline() {}
  virtual ~HorizontalVerticalBaseline() {}

  inline MapType operator ()(const vnl_vector<double> & baselineRTN)
  {
    MapType out;
	out.clear();
  	out.insert(std::pair<std::string,double>("Horizontal",baselineRTN(2)) );	
  	out.insert(std::pair<std::string,double>("Vertical",baselineRTN(1)) );	
	return out;
  }
};

} // end namespace otb::Functor


/** \class otbBaseline
 * Baseline is an abstract class for the baseline calculation
 *
 * \ingroup Operators
 */
template <class TMasterInputImage,class TSlaveInputImage, 
          class TBaselineFunctor= Functor::LengthOrientationBaseline>
class ITK_EXPORT Baseline : public itk::Object 
{
public:
  /** Standard class typedefs. */
  typedef Baseline                       Self;
  typedef itk::Object                        Superclass;
  typedef itk::SmartPointer<Self>            Pointer;
  typedef itk::SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Baseline, itk::Object);

  /** Type definition for the input image. */
  typedef TMasterInputImage  MasterImageType;
  typedef TSlaveInputImage   SlaveImageType;
  typedef TBaselineFunctor   FunctorType;

  /** Pointer type for the image. */
  typedef typename TMasterInputImage::Pointer  MasterImagePointer;
  typedef typename TSlaveInputImage::Pointer   SlaveImagePointer;
  
  /** Const Pointer type for the image. */
  typedef typename TMasterInputImage::ConstPointer MasterImageConstPointer;
  typedef typename TSlaveInputImage::ConstPointer SlaveImageConstPointer;

  typedef std::map<std::string,double> MapType;
  /** Set the input image. */
  void SetMasterPlateform(const ImageKeywordlist& image_kwl)
  {
	  	m_MasterPlateform->CreateSensorModel(image_kwl);
  }

  void SetSlavePlateform(const ImageKeywordlist& image_kwl)
  {
	  	m_SlavePlateform->CreateSensorModel(image_kwl); 
  }

  typedef otb::PlatformPositionAdapter     PlatformType;
  typedef typename PlatformType::Pointer   PlateformPointer;


  /** Compute the Baseline value. */
  virtual void Compute(double line);

  /** Get Master plateform position*/
  std::vector<double> GetMasterPlateformPosition(double line);

  /** Get Master plateform speed*/
  std::vector<double> GetMasterPlateformSpeed(double line);

  /** Get Slave plateform position*/
  std::vector<double> GetSlavePlateformPosition(double line);

  /** Get Slave plateform position*/
  std::vector<double> GetSlavePlateformSpeed(double line);

  /** Compute the Baseline value. */
  void EvaluateMasterAndSlavePosition(
	            double masterLine, double slaveLine,
				std::vector<double> & masterPosition,
				std::vector<double> & slavePosition);

  /** Compute the Baseline value. */
  void EvaluateMasterAndSlaveSpeed(
	            double masterLine, double slaveLine,
				std::vector<double> & masterSpeed,
				std::vector<double> & slaveSpeed);

  double GetBaselineValue(std::string name)
		{	
			MapType::const_iterator it;

			it = m_Baseline.find(name);
			double value = it->second; 
			return value;
		}

  MapType GetBaseline()
  {
	  return m_Baseline;
  }

vnl_vector<double> BaselineInRTNSystem(
				std::vector<double> & masterPosition,
				std::vector<double> & slavePosition,
				std::vector<double> & masterSpeed);

protected:
  Baseline();
  virtual ~Baseline() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  Baseline(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MapType                      m_Baseline;
  FunctorType                  m_Functor;
  PlateformPointer             m_MasterPlateform;
  PlateformPointer             m_SlavePlateform;
};

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbBaseline.txx"
#endif

#endif /* __otbBaseline_h */
