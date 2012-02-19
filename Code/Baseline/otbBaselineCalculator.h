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
#include <vnl/vnl_vector.h>

namespace otb
{

/** \class otbBaselineCalculator
 * Baseline is an abstract class for the baseline calculation
 *
 * \ingroup Operators
 */
template <class TMasterInputImage,class TSlaveInputImage, class TFunctor>
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

  /** Type definition for the input image. */
  typedef TMasterInputImage  MasterImageType;
  typedef TSlaveInputImage   SlaveImageType;

  /** Pointer type for the image. */
  typedef typename TMasterInputImage::Pointer  MasterImagePointer;
  typedef typename TSlaveInputImage::Pointer   SlaveImagePointer;
  
  /** Const Pointer type for the image. */
  typedef typename TMasterInputImage::ConstPointer MasterImageConstPointer;
  typedef typename TSlaveInputImage::ConstPointer SlaveImageConstPointer;

  /** Set the input image. */
  itkSetConstObjectMacro(MasterImage,MasterImageType);
  itkSetConstObjectMacro(SlaveImage,SlaveImageType);

  typedef otb::PlatformPositionToBaselineCalculator<TFunctor> BaselineType;
  typedef typename BaselineType::Pointer                      BaselinePointer;
  typedef typename BaselineType::ConstPointer                 BaselineConstPointer;
  typedef typename BaselineType::BaselineFunctorOutputType    OutputBaselineType;
  typedef typename BaselineType::BaselineCalculusEnumType     BaselineCalculusEnumType;

  /** Typedef for Coefficient */
  typedef vnl_vector<double>            CoefficientType;        

  /** Compute the Baseline value. */
  void Compute(BaselineCalculusEnumType map);

  double EvaluateBaseline(double row,double col);

  itkGetObjectMacro(PlateformPositionToBaselineCalculator,BaselineType);

  void ExtractBaseline(	BaselineCalculusEnumType map,
						std::vector<typename MasterImageType::PointType> & pointImage,
						std::vector<double> & baselineImage);


protected:
  BaselineCalculator();
  virtual ~BaselineCalculator() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;


  vnl_vector<double> BaselineLinearSolve(
						std::vector<typename MasterImageType::PointType> & pointImage,
						std::vector<double> & baselineImage);

private:
  BaselineCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MasterImageConstPointer      m_MasterImage;
  SlaveImageConstPointer       m_SlaveImage;
  CoefficientType              m_BaselineCoefficient;
  BaselinePointer              m_PlateformPositionToBaselineCalculator;
};

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbBaselineCalculator.txx"
#endif

#endif /* __otbBaseline_h */
