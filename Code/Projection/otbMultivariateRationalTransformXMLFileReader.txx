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
#ifndef __otbMultivariateRationalTransformXMLFileReader_txx
#define __otbMultivariateRationalTransformXMLFileReader_txx

#include "otbMultivariateRationalTransformXMLFileReader.h"
#include "itkMacro.h"
#include "itksys/SystemTools.hxx"
#include "tinyxml.h"

namespace otb {


template < class TScalarType, unsigned int Dimension>
MultivariateRationalTransformXMLFileReader<TScalarType,Dimension>
::MultivariateRationalTransformXMLFileReader(): m_FileName("")
{
  m_MultivariateRationalTransform = MultivariateRationalTransformType::New();
}

template < class TScalarType, unsigned int Dimension>
void
MultivariateRationalTransformXMLFileReader<TScalarType,Dimension>
::Read()
{
  // initialize the output landmark in case of multiple execution
  //m_MultivariateRationalTransform = MultivariateRationalTransformType::New();

  // Check if the filename is not empty
  if(m_FileName.empty())
    itkExceptionMacro(<<"The XML output FileName is empty, please set the filename via the method SetFileName");

  // Check that the right extension is given : expected .xml */
  if (itksys::SystemTools::GetFilenameLastExtension(m_FileName) != ".xml")
    {
    itkExceptionMacro(<<itksys::SystemTools::GetFilenameLastExtension(m_FileName)
                      <<" is a wrong Extension FileName : Expected .xml");
    }


  // Open the xml file
  TiXmlDocument doc(m_FileName.c_str());

  if (!doc.LoadFile())
    {
    itkExceptionMacro(<<"Can't open file "<<m_FileName);
    }
  
  TiXmlHandle hDoc(&doc);
  TiXmlHandle root    = hDoc.FirstChildElement("RATIONAL_TRANSFORM_FILE");

  int numeratorDegree = 0;
  int denominatorDegree = 0;
  int spaceDimension = 0;

  if(root.FirstChildElement("SPACE_DIMENSION").ToElement() == NULL)
    {
    itkExceptionMacro(<<"Needs the SPACE_DIMENSION tag in XML file");
    }
  root.FirstChildElement("SPACE_DIMENSION").ToElement()->QueryIntAttribute("VALUE", &spaceDimension);

  if(root.FirstChildElement("NUMERATOR_DEGREE").ToElement() == NULL)
    {
    itkExceptionMacro(<<"Needs the NUMERATOR_DEGREE tag in XML file");
    }
  root.FirstChildElement("NUMERATOR_DEGREE").ToElement()->QueryIntAttribute("VALUE", &numeratorDegree);

  if(root.FirstChildElement("DENOMINATOR_DEGREE").ToElement() == NULL)
    {
    itkExceptionMacro(<<"Needs the DENOMINATOR_DEGREE tag in XML file");
    }
  root.FirstChildElement("DENOMINATOR_DEGREE").ToElement()->QueryIntAttribute("VALUE", &denominatorDegree);

  if( numeratorDegree <0 )
    {
    itkExceptionMacro(<< "Numerator Degree must be greater or equal to zero");
    }
  if( denominatorDegree <0 )
    {
    itkExceptionMacro(<< "Denominator Degree must be greater or equal to zero");
    }

  if( spaceDimension !=m_MultivariateRationalTransform->SpaceDimension )
    {
    itkExceptionMacro(<< "spaceDimension must have the same dimension with the MultivariateRationalTransformXMLFileReader() dimension class");
    }

  m_MultivariateRationalTransform->SetNumeratorDegree( numeratorDegree );
  m_MultivariateRationalTransform->SetDenominatorDegree( denominatorDegree );

  typename MultivariateRationalTransformType::ParametersType params(m_MultivariateRationalTransform->GetNumberOfParameters());
  params.Fill(0);

  unsigned int strip = m_MultivariateRationalTransform->GetNumeratorDegree() + m_MultivariateRationalTransform->GetDenominatorDegree() + 2;

  for( TiXmlElement* currentDimension = root.FirstChildElement("COEFFICIENT_LIST").FirstChildElement("DIMENSION").ToElement();
      currentDimension != NULL;
      currentDimension = currentDimension->NextSiblingElement() )
        {
        // Get the value for Degree and polynomial attributes
        int dim;
         currentDimension->QueryIntAttribute("INDEX", &dim);
         if( (dim<0) || (dim > m_MultivariateRationalTransform->SpaceDimension) )
           {
           itkExceptionMacro(<< "index value range must be 0 to "<< m_MultivariateRationalTransform->SpaceDimension );
           }
         std::cout << "dim : " << dim<< std::endl;
         for( TiXmlElement* currentNumeratorCoef = currentDimension->FirstChildElement("NUMERATOR_COEFFICIENT")->FirstChildElement("DEGREE")->ToElement();
             currentNumeratorCoef != NULL;
             currentNumeratorCoef = currentNumeratorCoef->NextSiblingElement() )
           {
             int degree;
             double value;
             currentNumeratorCoef->QueryIntAttribute("INDEX", &degree);
             currentNumeratorCoef->QueryDoubleAttribute("VALUE", &value);
             if( (degree<0) || (degree > m_MultivariateRationalTransform->GetNumeratorDegree()) )
               {
               itkExceptionMacro(<< "index value range must be 0 to "<< m_MultivariateRationalTransform->GetNumeratorDegree() );
               }
             params[strip*dim + degree] = value;
           }

         for( TiXmlElement* currentDenominatorCoef = currentDimension->FirstChildElement("DENOMINATOR_COEFFICIENT")->FirstChildElement("DEGREE")->ToElement();
             currentDenominatorCoef != NULL;
             currentDenominatorCoef = currentDenominatorCoef->NextSiblingElement() )
           {
             int degree;
             double value;
             currentDenominatorCoef->QueryIntAttribute("INDEX", &degree);
             currentDenominatorCoef->QueryDoubleAttribute("VALUE", &value);
             if( (degree<0) || (degree > m_MultivariateRationalTransform->GetDenominatorDegree()) )
               {
               itkExceptionMacro(<< "index value range must be 0 to "<< m_MultivariateRationalTransform->GetDenominatorDegree() );
               }
             params[strip*dim+m_MultivariateRationalTransform->GetNumeratorDegree()+1 + degree] = value;
           }
        }
  m_MultivariateRationalTransform->SetParameters(params);
}

template < class TScalarType, unsigned int Dimension>
void
MultivariateRationalTransformXMLFileReader<TScalarType,Dimension>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  // Call superclass implementation
  Superclass::PrintSelf(os, indent);
  os << indent << "\nMultivariateRationalTransform : " << std::endl;
  m_MultivariateRationalTransform->Print(os);

}

} // End namespace otb

#endif
