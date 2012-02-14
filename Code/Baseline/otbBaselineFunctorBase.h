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
#ifndef __otbBaselineFunctorBase_h
#define __otbBaselineFunctorBase_h

#include "itkLightObject.h"
#include <vnl/vnl_vector.h>
#include "itkObject.h"

namespace otb
{

namespace Functor {

class ITK_EXPORT BaselineFunctorBase : public itk::Object
{
public:

  typedef BaselineFunctorBase                           Self;
  typedef itk::Object                                   Superclass;
  typedef itk::SmartPointer<Self>                       Pointer;
  typedef itk::SmartPointer<const Self>                 ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BaselineFunctorBase, itk::Object);

  typedef vnl_vector<double> VectorType;

  typedef enum {Horizontal, Vertical, Parallel, Perpendicular, 
					Length, Orientation } BaselineCalculusEnumType;

  virtual double GetHorizontalBaseline() const;
  virtual double GetVerticalBaseline() const;
  virtual double GetParallelBaseline() const;
  virtual double GetPerpendicularBaseline() const;
  virtual double GetLengthBaseline() const;
  virtual double GetOrientationBaseline() const;

  double GetBaseline(BaselineCalculusEnumType);

  void SetRTNBaseline(VectorType & RTNBaseline)
  {
	m_RTNBaseline = RTNBaseline;
  }

protected:
  BaselineFunctorBase() {};
  ~BaselineFunctorBase() {};

  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  BaselineFunctorBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  VectorType m_RTNBaseline;
};
	
} // end namespace functor

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbBaselineFunctorBase.cxx"
#endif


#endif /* __otbBaselineFunctorBase_h */
