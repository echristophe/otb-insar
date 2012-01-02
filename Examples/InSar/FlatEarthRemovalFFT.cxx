/*=========================================================================

   Copyright 2011 David Dubois
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

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbStreamingImageFileWriter.h"
#include "otbImageFileWriter.h"

#include "otbFlatEarthRemovalImageFilter.h"

// Command line:

// ./FlatEarthremovalFFT interf.tif remove_earth_interf.tif blockSize padSize

int main(int argc, char* argv[])
{

  if (argc != 5)
    {
    std::cerr << "Usage: " << argv[0] << " interferogram flatEarthRemoveInterferogram blockSize padSize" << std::endl;
    return EXIT_FAILURE;
    }

  typedef std::complex<double>    PixelType;
  typedef otb::Image<PixelType,2> ImageType;
  typedef double				  OutPixelType;
  typedef otb::Image<OutPixelType, 2> OutImageType;
  typedef ImageType::IndexType    IndexType;
  typedef double                  ScalarPixelType;
  typedef otb::Image<ScalarPixelType,2>                      ScalarImageType;
  
  typedef otb::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer interferogram = ReaderType::New();
  interferogram->SetFileName(argv[1]);
  interferogram->UpdateOutputInformation();

  typedef otb::FlatEarthRemovalImageFilter< ImageType, OutImageType > FlatEarthRemoveType;
  FlatEarthRemoveType::Pointer flatEarthRemove = FlatEarthRemoveType::New();
  flatEarthRemove->SetInput(interferogram->GetOutput());
  flatEarthRemove->SetPatchSizePerDim(atoi(argv[3]));
  flatEarthRemove->SetPadSizePerDim(atoi(argv[4]));

  typedef otb::StreamingImageFileWriter<OutImageType> WriterType;
  WriterType::Pointer flatEarthRemoveWriter = WriterType::New();
  flatEarthRemoveWriter->SetFileName(argv[2]);
  flatEarthRemoveWriter->SetInput(flatEarthRemove->GetOutput());
  flatEarthRemoveWriter->Update();


  return EXIT_SUCCESS;
}
