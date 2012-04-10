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
#ifndef __otbMultivariateRationalTransformXMLFileReader_h
#define __otbMultivariateRationalTransformXMLFileReader_h

#include "itkObject.h"
#include "otbObjectList.h"
#include "itkMacro.h"
#include "otbMultivariateRationalTransform.h"

namespace otb {

/** \class MultivariateRationalTransformXMLFileReader
 *  \brief   Read a xml file where are stored Deformation Model
 *           and build an MultivariateRationalTransform<TScalarType,2>

 *  Parameters in the parameters vector are in the following order:
 *  dim0num0 dim0num1 ... dim0numN dim0denom0 dim0denom1
 *  ... dim0denomM ... dim1num0 ... dimDdenomM.
 *
 * To get the output MultivariateRationalTransform<MultivariateRationalTransform<TScalarType,2>
 * use the method GetOutput()
 */
template < class TScalarType, unsigned int Dimension = 2>
class  MultivariateRationalTransformXMLFileReader :
    public itk::Object
{
public:
  /** Standard class typedefs */
  typedef MultivariateRationalTransformXMLFileReader    Self;
  typedef itk::Object                      Superclass;
  typedef itk::SmartPointer< Self >        Pointer;
  typedef itk::SmartPointer<const Self>    ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MultivariateRationalTransformXMLFileReader, itk::Object);
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

 /** PointSet convenient typedef */
  typedef TScalarType                                   ScalarType;
  

  /** input MultivariateRationalTransform type */
  typedef MultivariateRationalTransform<ScalarType, Dimension >    MultivariateRationalTransformType;
  typedef typename MultivariateRationalTransformType::Pointer      MultivariateRationalTransformPointer;
  
  /** Set the output filename */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  /** Get the output Landmark*/
  MultivariateRationalTransformPointer GetOutput()
  {
    return m_MultivariateRationalTransform;
  }

  /** Update*/
  void Update()
  {
    this->Read();
  }
  
protected:
  virtual void Read();

  MultivariateRationalTransformXMLFileReader();
  virtual ~MultivariateRationalTransformXMLFileReader() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  MultivariateRationalTransformXMLFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string                 m_FileName;
  MultivariateRationalTransformPointer    m_MultivariateRationalTransform;

}; // end of class MultivariateRationalTransformXMLFileReader

} // end of namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbMultivariateRationalTransformXMLFileReader.txx"
#endif

#endif
