/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See OTBCopyright.txt for details.


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __otbComplexInterpolateImageFunction_h
#define __otbComplexInterpolateImageFunction_h

#include "itkInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNumericTraits.h"

namespace otb
{

/** \class ComplexInterpolateImageFunction
 * \brief Generic interpolation of an otb::Image.
 *
 * ComplexInterpolateImageFunction interpolates image intensity according to a
 * resampling profil.
 *
 * The Initialize() method need to be call to create the filter.
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

  /** Dimension underlying input image. */
  //itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

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
  //unsigned int GetRadius() { return this->GetFunction().GetRadius(); };

  /** Set/Get the window radius*/
  // Don't have to be used here, just declared for the inheritance classes.
  //virtual void SetWindowSize(unsigned int win){ m_WindowSize = win; };

  /** Get the functor list */
  virtual FunctionType& GetFunction(void)
  {
    return m_Function;
  }

  /** Initialize tables: need to be call explicitely */
  virtual void Initialize();

  /** Weights normalization accessors*/
  itkSetMacro(NormalizeWeight, bool);
  itkGetMacro(NormalizeWeight, bool);

protected:
  ComplexInterpolateImageFunction();
  virtual ~ComplexInterpolateImageFunction();
  virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

  /** Call the superclass implementation and set the TablesHaveBeenGenerated
   * flag to false */
  virtual void Modified(void);

  /** Delete tables.*/
  virtual void ResetOffsetTable();
  /** Initialize used tables*/
  virtual void InitializeTables();
  /** Fill the weight offset table*/
  virtual void FillWeightOffsetTable();

private:
  ComplexInterpolateImageFunction(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented
  /** Store the window radius. */
  //unsigned int m_Radius;
  // Constant to store twice the radius
  unsigned int m_WindowSize;

  /** Used function */
  FunctionType m_Function;
  /** Store the image dimension.*/
  unsigned int m_ImageDimension;

  /** These members are declared mutable so that they can be
  regenerated seamlessly inside the EvaluateAtContinuousIndex method if
  they need to */
  /** Size of the offset table */
  mutable unsigned int m_OffsetTableSize;
  /** The offset array, used to keep a list of relevant
   * offsets in the neihborhoodIterator */
  mutable unsigned int *m_OffsetTable;
  /** Index into the weights array for each offset */
  mutable unsigned int **m_WeightOffsetTable;
  /** True if internal statistics have been generated */
  mutable bool m_TablesHaveBeenGenerated;
  /** Weights normalization */
  bool m_NormalizeWeight;

  double m_ZeroDoppler;
};

} // end namespace itk

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbComplexInterpolateImageFunction.txx"
#endif

#endif
