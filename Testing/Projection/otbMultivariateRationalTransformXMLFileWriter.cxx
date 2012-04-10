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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "otbMultivariateRationalTransformXMLFileWriter.h"
#include "otbMultivariateRationalTransform.h"


// typedef 

int otbMultivariateRationalTransformXMLFileWriter(int argc, char* argv[])
{
  typedef otb::MultivariateRationalTransformXMLFileWriter<double>  XMLFileWriterType;

  const char * outfname = argv[1];

  typedef otb::MultivariateRationalTransform<double> MultivariateRationalTransformType;
  MultivariateRationalTransformType::Pointer rt = MultivariateRationalTransformType::New();
  rt->SetNumeratorDegree(2);
  rt->SetDenominatorDegree(1);
  // Rational is
  // fx(x, y) = (1+2*x+3*x^2+4*x^3+5*x^4)/(6+7*x+8*x^2+9*x^3+10*x^4)
  // fy(x, y) = (11+12*y+13*y^2+14*y^3+15*y^4)/(16+17*y+18*y^2+19*y^3+20*y^4)

  MultivariateRationalTransformType::ParametersType params(rt->GetNumberOfParameters());
  // Rational is
  // fx(x, y) = (1+2y+3*y^2+(4+5y+6*y^2)*x + (7+8y+9*y^2)*x^2 )/(10+11*y+(12+13*y)*x)
  // fy(x, y) = (14+15y+16*y^2+(17+18y+19*y^2)*x + (20+21y+22*y^2)*x^2 )/(23+24*y+(25+26*y)*x)
  params[0]=1;
  params[1]=2;
  params[2]=3;
  params[3]=4;
  params[4]=5;
  params[5]=6;
  params[6]=7;
  params[7]=8;
  params[8]=9;
  params[9]=10;
  params[10]=11;
  params[11]=12;
  params[12]=13;
  params[13]=14;
  params[14]=15;
  params[15]=16;
  params[16]=17;
  params[17]=18;
  params[18]=19;
  params[19]=20;
  params[20]=21;
  params[21]=22;
  params[22]=23;
  params[23]=24;
  params[24]=25;
  params[25]=26;

  rt->SetParameters(params);

  // instanciation
  XMLFileWriterType::Pointer     xmlwriter = XMLFileWriterType::New();
  xmlwriter->SetFileName(outfname);
  xmlwriter->SetInput(rt);
  xmlwriter->Update();

  return EXIT_SUCCESS;
}
