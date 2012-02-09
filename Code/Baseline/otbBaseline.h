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

  inline MapType operator ()(const std::vector<double> & masterPosition,
							 const std::vector<double> & slavePosition	)
  {
    MapType out;
	out.clear();

	double baselineLength = vcl_sqrt(
					(masterPosition[0] - slavePosition[0]) * (masterPosition[0] - slavePosition[0]) +
					(masterPosition[1] - slavePosition[1]) * (masterPosition[1] - slavePosition[1]) +
					(masterPosition[2] - slavePosition[2]) * (masterPosition[2] - slavePosition[2]));
    //TODO 
  	out.insert(std::pair<std::string,double>("Length",baselineLength) );	
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

  inline MapType operator ()(const std::vector<double> & masterPosition,
							 const std::vector<double> & slavePosition	)
  {
    MapType out;
	out.clear();
	/*
	double baselineLength = vcl_sqrt(
					(masterPosition[0] - slavePosition[0]) * (masterPosition[0] - slavePosition[0]) +
					(masterPosition[1] - slavePosition[1]) * (masterPosition[1] - slavePosition[1]) +
					(masterPosition[2] - slavePosition[2]) * (masterPosition[2] - slavePosition[2]));

  	out.insert(std::pair<std::string,double>("Baseline",baselineLength) );	
    */
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

  inline MapType operator ()(const std::vector<double> & masterPosition,
							 const std::vector<double> & slavePosition	)
  {
    MapType out;
	out.clear();
	/*
	double baselineLength = vcl_sqrt(
					(masterPosition[0] - slavePosition[0]) * (masterPosition[0] - slavePosition[0]) +
					(masterPosition[1] - slavePosition[1]) * (masterPosition[1] - slavePosition[1]) +
					(masterPosition[2] - slavePosition[2]) * (masterPosition[2] - slavePosition[2]));

  	out.insert(std::pair<std::string,double>("Baseline",baselineLength) );	
    */
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
  itkSetConstObjectMacro(MasterImage,MasterImageType);
  itkSetConstObjectMacro(SlaveImage,SlaveImageType);



  /** Compute the Baseline value. */
  virtual void Compute(double line);

  /** Get Master plateform position*/
  std::vector<double> GetMasterPlateformPosition(double line);

  /** Get Slave plateform position*/
  std::vector<double> GetSlavePlateformPosition(double line);

  /** Compute the Baseline value. */
  void EvaluateMasterAndSlavePosition(
	            double masterLine, double slaveLine,
				std::vector<double> & masterPosition,
				std::vector<double> & slavePosition);

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


protected:
  Baseline();
  virtual ~Baseline() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  Baseline(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MasterImageConstPointer      m_MasterImage;
  SlaveImageConstPointer       m_SlaveImage;
  MapType                      m_Baseline;
  FunctorType                  m_Functor;

};

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbBaseline.txx"
#endif

#endif /* __otbBaseline_h */
