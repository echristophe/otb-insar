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
#ifndef __otbBaselineCalculator_h
#define __otbBaselineCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "otbPlatformPositionToBaselineCalculator.h"
#include "itkImageRegion.h"
#include "itkPoint.h"
#include <vnl/vnl_vector.h>

namespace otb
{

/** \class otbBaselineCalculator
 * Baseline is an abstract class for the baseline calculation
 *
 * \ingroup Operators
 */
template <class TFunctor,unsigned int Dimension = 2>
class ITK_EXPORT BaselineCalculator : public itk::Object 
{
public:
  /** Standard class typedefs. */
  typedef BaselineCalculator                 Self;
  typedef itk::Object                        Superclass;
  typedef itk::SmartPointer<Self>            Pointer;
  typedef itk::SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BaselineCalculator, itk::Object);

  itkStaticConstMacro(SpaceDimension, unsigned int, Dimension);


  typedef otb::PlatformPositionToBaselineCalculator<TFunctor>                              PlateformPositionToBaselineCalculatorType;
  typedef typename PlateformPositionToBaselineCalculatorType::Pointer                      PlateformPositionToBaselinePointer;
  typedef typename PlateformPositionToBaselineCalculatorType::ConstPointer                 BaselineConstPointer;
  typedef typename PlateformPositionToBaselineCalculatorType::BaselineFunctorOutputType    OutputBaselineType;
  typedef typename PlateformPositionToBaselineCalculatorType::BaselineCalculusEnumType     BaselineCalculusEnumType;

  typedef itk::ImageRegion< Dimension >			ImageRegionType;

  typedef itk::Point< double, Dimension >       PointType; 

  /** Typedef for Coefficient */
  typedef vnl_vector<double>            CoefficientType;        

  /** Compute the Baseline value. */
  void Compute(BaselineCalculusEnumType map);

  double EvaluateBaseline(double row,double col);

  itkSetObjectMacro(PlateformPositionToBaselineCalculator,PlateformPositionToBaselineCalculatorType);
  itkGetObjectMacro(PlateformPositionToBaselineCalculator,PlateformPositionToBaselineCalculatorType);

  itkSetMacro(Region,ImageRegionType);
  itkGetConstMacro(Region,ImageRegionType);

  itkSetMacro(LineOffsetWithMaster,double);
  itkGetConstMacro(LineOffsetWithMaster,double);

  void ExtractBaseline(	BaselineCalculusEnumType map,
						std::vector<PointType> & pointImage,
						std::vector<double> & baselineImage);

  void FoundMinimumTangentialBaseline();

protected:
  BaselineCalculator();
  virtual ~BaselineCalculator() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;


  vnl_vector<double> BaselineLinearSolve(
						std::vector<PointType> & pointImage,
						std::vector<double> & baselineImage);

private:
  BaselineCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ImageRegionType                     m_Region;
  CoefficientType                     m_BaselineCoefficient;
  PlateformPositionToBaselinePointer  m_PlateformPositionToBaselineCalculator;
  double							  m_LineOffsetWithMaster;	
};

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbBaselineCalculator.txx"
#endif

#endif /* __otbBaseline_h */
