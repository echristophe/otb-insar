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


#include "otbMultivariateRationalTransform.h"
#include <fstream>

int otbMultivariateRationalTransformNew(int argc, char* argv[])
{
  typedef otb::MultivariateRationalTransform<> MultivariateRationalTransformType;

  // Instantiation
  MultivariateRationalTransformType::Pointer mrt = MultivariateRationalTransformType::New();

  return EXIT_SUCCESS;
}

int otbMultivariateRationalTransformTest(int argc, char* argv[])
{
  typedef otb::MultivariateRationalTransform<> MultivariateRationalTransformType;

  // Instantiation
  MultivariateRationalTransformType::Pointer mrt = MultivariateRationalTransformType::New();
  mrt->SetNumeratorDegree(2);
  mrt->SetDenominatorDegree(1);

  MultivariateRationalTransformType::ParametersType params(mrt->GetNumberOfParameters());
  params.Fill(1.);

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

  mrt->SetParameters(params);

  MultivariateRationalTransformType::InputPointType inputPoint;
  MultivariateRationalTransformType::OutputPointType outputPoint;

  std::ofstream ofs;
  ofs.open(argv[1]);

  // Set floatfield to format writing properly
  ofs.setf(std::ios::fixed, std::ios::floatfield);
  ofs.precision(10);

  unsigned int idx = 2;

  ofs<<"Multivariate Rational function is: "<<std::endl;
  ofs<< "fx(x, y) = (1+2y+3*y^2+(4+5y+6*y^2)*x + (7+8y+9*y^2)*x^2 )/(6+7*y+(8+9*y)*x)"<<std::endl;
  ofs<< "fy(x, y) = (10+11y+12*y^2+(13+14y+15*y^2)*x + (16+17y+18*y^2)*x^2 )/(19+20*y+(21+22*y)*x)"<<std::endl;

  while(idx+1<(unsigned int)argc)
    {
    inputPoint[0] = atof(argv[idx]);
    inputPoint[1] = atof(argv[idx+1]);
    outputPoint = mrt->TransformPoint(inputPoint);
    ofs<<inputPoint<<" -> "<<outputPoint<<std::endl;
    idx+=2;
    }

  ofs.close();

  return EXIT_SUCCESS;
}

