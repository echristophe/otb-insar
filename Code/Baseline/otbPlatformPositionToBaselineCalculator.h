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
#ifndef __otbPlatformPositionToBaselineCalculator_h
#define __otbPlatformPositionToBaselineCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "otbMacro.h"
#include <vnl/vnl_vector.h>
#include "otbPlatformPositionAdapter.h"
#include "otbBaselineFunctorBase.h"

namespace otb
{
/** \class otbPlatformPositionToBaselineCalculator
 * Baseline is an abstract class for the baseline calculation
 *
 * \ingroup Operators
 */
class ITK_EXPORT PlatformPositionToBaselineCalculator : public itk::Object 
{
public:
  /** Standard class typedefs. */
  typedef PlatformPositionToBaselineCalculator  Self;
  typedef itk::Object                           Superclass;
  typedef itk::SmartPointer<Self>               Pointer;
  typedef itk::SmartPointer<const Self>         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PlatformPositionToBaselineCalculator, itk::Object);

  /** Type definition for Plateform Position adapter. */
  typedef otb::PlatformPositionAdapter     PlatformType;
  typedef PlatformType::Pointer            PlateformPointer;

  /** Type definition for Baseline Functor. */
  typedef Functor::BaselineFunctorBase   BaselineFunctorType;
  typedef BaselineFunctorType::Pointer   BaselineFunctorPointer;


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


  /**
   * Set/Get the BaselineFunctor object.
   */
  BaselineFunctorPointer GetBaselineFunctor() { return m_BaselineFunctor; }

  void SetBaselineFunctor( BaselineFunctorType * baselineFunctor )
    {
    if ( m_BaselineFunctor != baselineFunctor )
      {
      m_BaselineFunctor = baselineFunctor;
      this->Modified();
      }
    }

  /**
   * Enum type that provides for an easy interface to existing BaselineFunctor.
   */
  typedef enum {	HorizontalVertical, ParallelPerpendicular, 
					LengthOrientation } BaselineFunctorEnumType;

  void SetBaselineFunctor( BaselineFunctorEnumType );


  /** Compute the Baseline value. */
  virtual void Compute(double line);

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
  PlatformPositionToBaselineCalculator();
  virtual ~PlatformPositionToBaselineCalculator() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  PlatformPositionToBaselineCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MapType                          m_Baseline;
  BaselineFunctorPointer           m_BaselineFunctor;
  PlateformPointer                 m_MasterPlateform;
  PlateformPointer                 m_SlavePlateform;
};

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbPlatformPositionToBaselineCalculator.cxx"
#endif

#endif /* __otbBaseline_h */
