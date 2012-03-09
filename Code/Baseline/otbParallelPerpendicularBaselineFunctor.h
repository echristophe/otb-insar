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
#ifndef __otbParallelPerpendicularBaselineFunctor_h
#define __otbParallelPerpendicularBaselineFunctor_h

#include "itkObjectFactory.h"
#include "itkLightObject.h"
#include "otbBaselineFunctorBase.h"
#include <vnl/vnl_vector.h>

namespace otb
{

namespace Functor {

class ITK_EXPORT ParallelPerpendicularBaselineFunctor : public BaselineFunctorBase
{
public:

  typedef ParallelPerpendicularBaselineFunctor          Self;
  typedef BaselineFunctorBase                           Superclass;
  typedef itk::SmartPointer<Self>                       Pointer;
  typedef itk::SmartPointer<const Self>                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro(ParallelPerpendicularBaselineFunctor, BaselineFunctorBase);

  typedef Superclass::VectorType     VectorType;
  typedef double                     RealType;

  typedef enum {Parallel, Perpendicular, Tangential} BaselineCalculusEnumType;

  OutputType GetBaseline(BaselineCalculusEnumType);

  virtual OutputType GetParallelBaseline() const;
  virtual OutputType GetPerpendicularBaseline() const;

  itkSetMacro( LookDirection, RealType );
  itkGetConstMacro( LookDirection, RealType  );

protected:
  ParallelPerpendicularBaselineFunctor() {};
  ~ParallelPerpendicularBaselineFunctor() {};


void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  ParallelPerpendicularBaselineFunctor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double m_LookDirection;

};

		
} // end namespace functor
} // end namespace otb

#endif /* __otbParallelPerpendicularBaselineFunctor_h */
