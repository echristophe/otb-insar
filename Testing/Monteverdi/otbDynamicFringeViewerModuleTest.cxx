/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See OTBCopyright.txt for details.


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "otbDynamicFringeViewerModule.h"

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"

int otbDynamicFringeViewerModuleTest(int argc, char* argv[])
{

  otb::DynamicFringeViewerModule::Pointer specificModule = otb::DynamicFringeViewerModule::New();
  otb::Module::Pointer module = specificModule.GetPointer();
  
  std::cout<<"Module: "<<module<<std::endl;

  // Put in the tests
  const char * infname = argv[1];
  bool  run = atoi(argv[2]);
  
  typedef otb::Image<double,2>                ImageType;
  typedef otb::ImageFileReader<ImageType>     ReaderType;
  typedef otb::ImageFileWriter<ImageType>     WriterType;
  
  // Reader
  ReaderType::Pointer reader = ReaderType::New();
  
  reader->SetFileName(infname);
  reader->GenerateOutputInformation();

  otb::DataObjectWrapper wrapperIn = otb::DataObjectWrapper::Create(reader->GetOutput());
#if 0
  std::cout<<"Input wrapper: "<<wrapperIn<<std::endl;
  
  module->AddInputByKey("FirstImage",wrapperIn);
  module->Start();
  
  // Simulate Addition Operation 
  specificModule->guiOperation->value(0);
  specificModule->bOk->do_callback();

  
  
  if(run)
    {
      Fl::run();
    }
  else
    {
      Fl::check();
    }


  otb::DataObjectWrapper wrapperOut = module->GetOutputByKey("OutputImage");
  std::cout<<"Output wrapper: "<<wrapperOut<<std::endl;
  ImageType::Pointer outImage = dynamic_cast<ImageType *>(wrapperOut.GetDataObject());
  
  //Write the image
  WriterType::Pointer  writer = WriterType::New();
  writer->SetFileName(argv[3]);
  writer->SetInput(outImage);
  writer->Update();
#endif
  return EXIT_SUCCESS;

}

