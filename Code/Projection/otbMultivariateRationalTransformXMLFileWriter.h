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
#ifndef __otbMultivariateRationalTransformXMLFileWriter_h
#define __otbMultivariateRationalTransformXMLFileWriter_h

#include "itkProcessObject.h"
#include "otbMultivariateRationalTransform.h"
#include <utility>
#include <string>

namespace otb {

/** \class MultivariateRationalTransformXMLFileWriter
 *  \brief Write in a xml file the deformation model
 *
 */
template < class TScalarType, unsigned int Dimension = 2>
class  MultivariateRationalTransformXMLFileWriter :
    public itk::Object
{
public:
  /** Standard class typedefs */
  typedef MultivariateRationalTransformXMLFileWriter    Self;
  typedef itk::Object                      Superclass;
  typedef itk::SmartPointer< Self >        Pointer;
  typedef itk::SmartPointer<const Self>    ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MultivariateRationalTransformXMLFileWriter, itk::Object);
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** ScalarType convenient typedef */
   typedef TScalarType                                   ScalarType;


   /** input MultivariateRationalTransform type */
   typedef MultivariateRationalTransform<ScalarType, Dimension >    MultivariateRationalTransformType;
   typedef typename MultivariateRationalTransformType::Pointer      MultivariateRationalTransformPointer;
   typedef typename MultivariateRationalTransformType::ParametersType      MultivariateRationalTransformParametersType;

  /** Method to set/get the reference PointSet*/
  void SetInput(const MultivariateRationalTransformPointer MultivariateRationalTransform );
//  void SetInput(const MultivariateRationalTransformParametersType & MultivariateRationalTransform );

  /** Trigger the processing */
  void Update()
  {
    this->GenerateData();
  }
  
  /** Set the output filename */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);
  
protected:

  virtual void GenerateData();

  MultivariateRationalTransformXMLFileWriter();
  virtual ~MultivariateRationalTransformXMLFileWriter() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  MultivariateRationalTransformXMLFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string                m_FileName;
  MultivariateRationalTransformPointer   m_MultivariateRationalTransform;
}; // end of class MultivariateRationalTransformXMLFileWriter

} // end of namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbMultivariateRationalTransformXMLFileWriter.txx"
#endif

#endif
