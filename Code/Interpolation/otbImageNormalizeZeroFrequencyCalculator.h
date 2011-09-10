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

#ifndef __otbImageNormalizeZeroFrequencyCalculator_h
#define __otbImageNormalizeZeroFrequencyCalculator_h

#include "otbMacro.h"
#include "otbImage.h"
#include "itkNumericTraits.h"

namespace otb
{

/** \class ImageNormalizeZeroFrequencyCalculator
 * \brief Compute Normalize zero frequency.
 *
 * This class provides methods for computing the normalize zero frequency
 * of a single image.  
 *
 *
 * \ingroup Operators
 *
 */
template < class TImage >
class ITK_EXPORT ImageNormalizeZeroFrequencyCalculator : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ImageNormalizeZeroFrequencyCalculator<TImage>   Self;
  typedef itk::Object                                     Superclass;
  typedef itk::SmartPointer<Self>                         Pointer;
  typedef itk::SmartPointer<const Self>                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageNormalizeZeroFrequencyCalculator, itk::Object);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImage::ImageDimension);

  /** Standard scalar type within this class. */
  typedef double                       ScalarType;

  /** Standard vector type within this class. */
  typedef itk::Vector<ScalarType,itkGetStaticConstMacro(ImageDimension)> VectorType;

  /** Standard image type within this class. */
  typedef TImage ImageType;

  /** RealType typedef support. */
  typedef typename itk::NumericTraits<typename ImageType::PixelType>::ValueType ValueType;

  /** Standard image type pointer within this class. */
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename ImageType::PixelType    PixelType;

  typedef double  ScalarPixelType;
  typedef otb::Image<ScalarPixelType,ImageDimension> ScalarImageType;

  /** Set the input image. */
  virtual void SetImage( const ImageType * image )
    {
    if ( m_Image != image )
      {
      m_Image = image;
      this->Modified();
      m_Valid = false;
      }
    }



  /** Compute normalize zero frequency of a new or modified image.
   * This method computes the normalize zero frequency of the image given 
   */
  void Compute( void );

  /** Return the Normalize Zero Frequency. */
  VectorType GetNormalizeZeroFrequency() const;

protected:
  ImageNormalizeZeroFrequencyCalculator();
  virtual ~ImageNormalizeZeroFrequencyCalculator();
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  ImageNormalizeZeroFrequencyCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool       m_Valid;                     // Have Normalize Zero Frequency been computed yet?
  VectorType m_NormalizeZeroFrequency;    // Normalize Zero Frequency

  ImageConstPointer         m_Image;

};  // class ImageNormalizeZeroFrequencyCalculator

} // end namespace otb


#ifndef ITK_MANUAL_INSTANTIATION
#include "otbImageNormalizeZeroFrequencyCalculator.txx"
#endif

#endif /* __otbImageNormalizeZeroFrequencyCalculator_h */
