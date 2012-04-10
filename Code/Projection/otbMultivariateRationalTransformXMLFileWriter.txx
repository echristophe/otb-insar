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
#ifndef __otbMultivariateRationalTransformXMLFileWriter_txx
#define __otbMultivariateRationalTransformXMLFileWriter_txx

#include "otbMultivariateRationalTransformXMLFileWriter.h"
#include "itkMacro.h"
#include "itksys/SystemTools.hxx"
#include "tinyxml.h"

namespace otb {

template < class TScalarType, unsigned int Dimension>
MultivariateRationalTransformXMLFileWriter<TScalarType,Dimension>
::MultivariateRationalTransformXMLFileWriter(): m_FileName("")
{}


template < class TScalarType, unsigned int Dimension>
void
MultivariateRationalTransformXMLFileWriter<TScalarType,Dimension>
::SetInput(const MultivariateRationalTransformPointer inputMultivariateRationalTransform )
{
  m_MultivariateRationalTransform = inputMultivariateRationalTransform;
}

//void SetInput(const MultivariateRationalTransformParametersType & MultivariateRationalTransform );


template < class TScalarType, unsigned int Dimension>
void
MultivariateRationalTransformXMLFileWriter<TScalarType,Dimension>
::GenerateData()
{
  // Check if the filename is not empty
  if(m_FileName.empty())
    itkExceptionMacro(<<"The XML output FileName is empty, please set the filename via the method SetFileName");
  
  // Check that the right extension is given : expected .xml */
  if (itksys::SystemTools::GetFilenameLastExtension(m_FileName) != ".xml")
    {
    itkExceptionMacro(<<itksys::SystemTools::GetFilenameLastExtension(m_FileName)
                      <<" is a wrong Extension FileName : Expected .xml");
    }
  
  // Write the XML file
  TiXmlDocument doc;

  TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
  doc.LinkEndChild( decl );

  TiXmlElement * root = new TiXmlElement( "RATIONAL_TRANSFORM_FILE");
  doc.LinkEndChild( root );


  // Create a feature NUMERATOR_DEGREE
  TiXmlElement * spaceDimension = new TiXmlElement("SPACE_DIMENSION");
  spaceDimension->SetDoubleAttribute("VALUE", m_MultivariateRationalTransform->SpaceDimension);
  root->LinkEndChild( spaceDimension );

  // Create a feature NUMERATOR_DEGREE
  TiXmlElement * numeratorDegree = new TiXmlElement("NUMERATOR_DEGREE");
  numeratorDegree->SetDoubleAttribute("VALUE", m_MultivariateRationalTransform->GetNumeratorDegree());
  root->LinkEndChild( numeratorDegree );

  // Create a feature DENOMINATOR_DEGREE
  TiXmlElement * denominatorDegree = new TiXmlElement("DENOMINATOR_DEGREE");
  denominatorDegree->SetDoubleAttribute("VALUE", m_MultivariateRationalTransform->GetDenominatorDegree());
  root->LinkEndChild( denominatorDegree );


  // Create a feature PARAMETER_LIST
  TiXmlElement * parameterList = new TiXmlElement("COEFFICIENT_LIST");
  root->LinkEndChild( parameterList );

  typename MultivariateRationalTransformType::ParametersType params(m_MultivariateRationalTransform->GetNumberOfParameters());
  params = m_MultivariateRationalTransform->GetParameters();

  unsigned int strip = m_MultivariateRationalTransform->GetNumeratorDegree() + m_MultivariateRationalTransform->GetDenominatorDegree() + 2;
  for(unsigned int dim = 0 ; dim < m_MultivariateRationalTransform->SpaceDimension ; dim++)
    {
    TiXmlElement * dimension = new TiXmlElement("DIMENSION");
    dimension->SetAttribute("INDEX", dim);
    parameterList->LinkEndChild(dimension);
    if(m_MultivariateRationalTransform->GetNumeratorDegree()>0)
      {
          TiXmlElement * parameterNumerator = new TiXmlElement("NUMERATOR_COEFFICIENT");
          dimension->LinkEndChild( parameterNumerator );
          for(unsigned int index = 0 ; index <= m_MultivariateRationalTransform->GetNumeratorDegree() ; index++)
            {
            TiXmlElement * parameter = new TiXmlElement("DEGREE");
            parameter->SetAttribute("INDEX", index);
            parameter->SetDoubleAttribute("VALUE", params[strip*dim+index]);
            parameterNumerator->LinkEndChild(parameter);
            }
      }

    if(m_MultivariateRationalTransform->GetDenominatorDegree()>0)
      {
            TiXmlElement * parameterDenominator = new TiXmlElement("DENOMINATOR_COEFFICIENT");
            dimension->LinkEndChild( parameterDenominator );
            for(unsigned int index = 0 ; index <= m_MultivariateRationalTransform->GetNumeratorDegree() ; index++)
              {
              TiXmlElement * parameter = new TiXmlElement("DEGREE");
              parameter->SetAttribute("INDEX", index);
              parameter->SetDoubleAttribute("VALUE", params[strip*dim+m_MultivariateRationalTransform->GetNumeratorDegree()+1 +index]);
              parameterDenominator->LinkEndChild(parameter);
              }
      }
    }

  // Finally, write the file
  doc.SaveFile( m_FileName.c_str() );
}


template < class TScalarType, unsigned int Dimension>
void
MultivariateRationalTransformXMLFileWriter<TScalarType,Dimension>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  // Call superclass implementation
  Superclass::PrintSelf(os, indent);
  
}

} // End namespace otb

#endif
