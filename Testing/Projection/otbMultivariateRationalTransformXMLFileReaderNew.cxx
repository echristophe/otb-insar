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

#include "otbMultivariateRationalTransformXMLFileReader.h"

int otbMultivariateRationalTransformXMLFileReaderNew(int argc, char* argv[])
{
  typedef otb::MultivariateRationalTransformXMLFileReader<double>  MultivariateRationalTransformXMLFileReaderType;
  MultivariateRationalTransformXMLFileReaderType::Pointer     xmlreader = MultivariateRationalTransformXMLFileReaderType::New();
  return EXIT_SUCCESS;
}
