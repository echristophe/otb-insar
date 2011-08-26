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
#ifndef __otbComplexInterpolateImageFunction_h
#define __otbComplexInterpolateImageFunction_h

#include "itkInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNumericTraits.h"
#include "itkVariableLengthVector.h"

namespace otb
{

/** \class ComplexInterpolateImageFunction
 * \brief Complex Generic interpolation of an otb::Image.
 *
 * ComplexInterpolateImageFunction interpolates complex image according to a
 * resampling profil.
 *
 * Along the resampling profil, the processing steps are  :  
 *		- Evaluating Kernel value at the given position                    
 * 		- Evaluating the bandshifted phase term: the frequency offset from baseband. 
 *		- Apply the bandshifted phase term to the PixelValue
 * 
 * And with the previous result applies a back bandshiftedifted phase term. 
 * Hence, value is still frequency offseted from the baseband.
 *
 *
 * \ingroup ImageFunctions ImageInterpolators
 */

template <class TInputImage, class TFunction, class TBoundaryCondition = itk::ZeroFluxNeumannBoundaryCondition<TInputImage>,
    class TCoordRep = double>
class ITK_EXPORT ComplexInterpolateImageFunction :
  public itk::InterpolateImageFunction<TInputImage, TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef ComplexInterpolateImageFunction                       Self;
  typedef itk::InterpolateImageFunction<TInputImage, TCoordRep> Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ComplexInterpolateImageFunction, itk::InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Input and output images typedef definition. */
  typedef typename Superclass::OutputType     OutputType;
  typedef typename Superclass::InputImageType InputImageType;

  /** Index and typedef support. */
  typedef typename Superclass::IndexType                                     IndexType;
  typedef typename InputImageType::SizeType                                  SizeType;
  typedef TFunction                                                          FunctionType;
  typedef itk::ConstNeighborhoodIterator<InputImageType, TBoundaryCondition> IteratorType;

  /** RealType typedef support. */
  typedef typename itk::NumericTraits<typename InputImageType::PixelType>::RealType RealType;
  
  /** ScalarRealType typedef support. */
  typedef typename itk::NumericTraits<typename InputImageType::PixelType>::ScalarRealType ScalarRealType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** VariableLengthVectorType typedef support*/
  typedef itk::VariableLengthVector<double> VariableLengthVectorType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType& index) const;

  /** Set/Get the window radius*/
  virtual void SetRadius(unsigned int rad);
  virtual unsigned int GetRadius() const
  {
    return m_Function.GetRadius();
  }

  /** Get the functor list */
  virtual FunctionType& GetFunction(void)
  {
    return m_Function;
  }

  /** Weights normalization accessors*/
  itkSetMacro(NormalizeWeight, bool);
  itkGetConstMacro(NormalizeWeight, bool);

  /** Sum Normalize Weight  accessors*/
  itkSetMacro(IsSumNormalizeWeight, bool);
  itkGetConstMacro(IsSumNormalizeWeight, bool);
  
  /** Normalize Frequency offset from baseband accessors */
  itkGetConstMacro(NormalizeZeroFrequency,VariableLengthVectorType);

  void SetNormalizeZeroFrequency(double value)
  {
	m_NormalizeZeroFrequency.Fill(value);
  }

  void SetNormalizeZeroFrequency(VariableLengthVectorType value)
  {
	  if(value.Size() != m_NormalizeZeroFrequency.Size())
	  {
		itkExceptionMacro(<< "The agrument dimension of the SetNormalizeZeroFrequency() method must be the same dimension as the input image");
	  }
	  m_NormalizeZeroFrequency = value ;
  }


protected:
  ComplexInterpolateImageFunction();
  virtual ~ComplexInterpolateImageFunction();
  virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  ComplexInterpolateImageFunction(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  /** Used function */
  FunctionType m_Function;

  /** Weights normalization */
  bool m_NormalizeWeight;

  /** Normalization with the sum (L1 norm) or square (L2 norm) of the resampling profil*/
  bool m_IsSumNormalizeWeight;

  /** Normalize Frequency offset from baseband*/
  VariableLengthVectorType m_NormalizeZeroFrequency;
};

} // end namespace itk

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbComplexInterpolateImageFunction.txx"
#endif

#endif
