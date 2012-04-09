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
#ifndef __otbModuloImageFilter_h
#define __otbModuloImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace otb
{
  
/** \class ModuloImageFilter
 * \brief Computes the function Modulo pixel-wise
 *
 * Every output pixel is equal to (x+a) Modulo K . where x is the intensity of the
 * homologous input pixel, a us a offset provided by user and K is a user-provided constant.
 * 
 * \ingroup IntensityImageFilters  Multithreaded
 *
 */

namespace Function {  
template< class TInput, class TOutput>
class Modulo
{
public:
  Modulo() { m_Factor = 1.0; m_Offset =0.0; }
  ~Modulo() {};

  bool operator!=( const Modulo & other ) const
    {
		if( (m_Factor != other.m_Factor) && (m_Offset != other.m_Offset) )
      {
      return true;
      }
    return false;
    }
  bool operator==( const Modulo & other ) const
    {
    return !(*this != other);
    }
  
  inline TOutput operator()( const TInput & A ) const
    {
    double tempValue = static_cast<double>(A) - m_Offset;
    return static_cast<TOutput>( tempValue - floor(static_cast<double>(tempValue) /m_Factor)* m_Factor );
    }

  void SetFactor( double factor )
    {
    m_Factor = factor;
    }
  double GetFactor() const
    {
    return m_Factor;
    }
  void SetOffset( double offset )
    {
    m_Offset = offset;
    }
  double GetOffset() const
    {
    return m_Offset;
    }
private:
  double  m_Factor;
  double  m_Offset;
}; 
}
template <class TInputImage, class TOutputImage>
class ITK_EXPORT ModuloImageFilter :
    public
itk::UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        otb::Function::Modulo< 
  typename TInputImage::PixelType, 
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef ModuloImageFilter                                Self;
  typedef itk::UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    otb::Function::Modulo< typename TInputImage::PixelType, 
                           typename TOutputImage::PixelType> >  Superclass;
  typedef itk::SmartPointer<Self>                                    Pointer;
  typedef itk::SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Runtime information support. */
  itkTypeMacro(ModuloImageFilter, 
               UnaryFunctorImageFilter);

  void SetFactor( double factor )
    {
    if( factor == this->GetFunctor().GetFactor() ) 
      {
      return;
      }
    this->GetFunctor().SetFactor( factor );
    this->Modified();
    }


  void SetOffset( double offset )
    {
    if( offset == this->GetFunctor().GetOffset() ) 
      {
      return;
      }
    this->GetFunctor().SetOffset( offset );
    this->Modified();
    }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToDoubleCheck,
    (itk::Concept::Convertible<typename TInputImage::PixelType, double>));
  itkConceptMacro(DoubleConvertibleToOutputCheck,
    (itk::Concept::Convertible<double, typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  ModuloImageFilter() {}
  virtual ~ModuloImageFilter() {}

private:
  ModuloImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace otb


#endif
