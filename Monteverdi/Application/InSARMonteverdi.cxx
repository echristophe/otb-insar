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

#include "ConfigureMonteverdi.h"

#include "otbMonteverdiModel.h"
#include "otbMonteverdiViewGUI.h"
#include "otbMonteverdiController.h"
#include "otbSplashScreen.h"
#include "otbI18n.h"
#include "otbMsgReporter.h"

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_PNG_Image.H>

#include <string>
#include <ctime>
#include <iostream>


// There are function prototype conflits under cygwin between standard w32 API
// and standard C ones
#ifndef CALLBACK
#if defined(_WINDOWS) || defined(__CYGWIN__)
#define CALLBACK __stdcall
#else
#define CALLBACK
#endif
#endif

#include "otbReaderModule.h"
#include "otbWriterModule.h"
#include "otbWriterMVCModule.h"
#include "otbExtractROIModule.h"
#include "otbViewerModule.h"
#include "otbCachingModule.h"

//InSAR Modules
#include "otbDynamicFringeViewerModule.h"

int main(int argc, char* argv[])
{
  //Internationalization
  otbI18nMacro();


  // Parse command line parameters
  typedef otb::CommandLineArgumentParser ParserType;
  ParserType::Pointer parser = ParserType::New();
  parser->AddOption("--InputImage", "input image file name to be visualized in the monteverdi viewer", "-in", 1, false);
  parser->AddOptionNParams("--InputList", "inputs can be images and vectorDatas (are not visualized in the viewer)", "-il", false);
  parser->SetProgramDescription("InSarMonteverdi launcher");
  //   parser->AddOption("--NoSplashScreen", "Deactivate the splash screen","-NoSplash", 0, false);

  typedef otb::CommandLineArgumentParseResult ParserResultType;
  ParserResultType::Pointer parseResult = ParserResultType::New();

  try
    {
    parser->ParseCommandLine(argc, argv, parseResult);
    }
  catch (itk::ExceptionObject& err)
    {
    std::string descriptionException = err.GetDescription();
    if (descriptionException.find("ParseCommandLine(): Help Parser") != std::string::npos)
      {
      return EXIT_SUCCESS;
      }
    }

  // Splash Screen
  typedef otb::SplashScreen::Pointer SplashScreenPointerType;

 // SplashScreenPointerType splash = otb::SplashScreen::New();
 // splash->Build();
 // splash->Show();


  // Application
  typedef otb::MonteverdiModel       ModelType;
  typedef otb::MonteverdiController  ControllerType;
  typedef otb::MonteverdiViewGUI     ViewType;

  // Create the MVC
  ModelType::Pointer model = otb::MonteverdiModel::GetInstance();
  ViewType::Pointer view = ViewType::New();
  ControllerType::Pointer controller = ControllerType::New();
  controller->SetView(view);
  view->SetMonteverdiController(controller);

  // Register modules
  model->RegisterModule<otb::ReaderModule>("Reader","File/Open dataset");
  model->RegisterModule<otb::WriterModule> ("Writer","File/Save dataset");
  model->RegisterModule<otb::ExtractROIModule>("ExtractROI","File/Extract ROI from dataset");
  model->RegisterModule<otb::WriterMVCModule> ("Specific writer for X image","File/Save dataset (advanced)");
  model->RegisterModule<otb::CachingModule>("zCaching","File/Cache dataset");
//  model->RegisterModule<otb::ConcatenateModule>("Concatenate","File/Concatenate images");

  model->RegisterModule<otb::ViewerModule>("Viewer","Visualization/Viewer");


  // InSAR Module to add in monteverdi
  model->RegisterModule<otb::DynamicFringeViewerModule>("Dynamic Fringe Viewer","Visualization/InSAR Fringe Viewer");
//  model->RegisterModule<otb::ObjectLabelingOSMModule>("OSM Object Labeling","Learning/OSM OBIA");
  
//  model->RegisterModule<otb::MeanShiftModule> ("MeanShift","Filtering/Mean shift clustering");
//  model->RegisterModule<otb::PanSharpeningModule> ("PanSharpening","Filtering/Pan-sharpen an image");
//  model->RegisterModule<otb::FeatureExtractionModule>("FeatureExtraction", "Filtering/Feature extraction");
  
//  model->RegisterModule<otb::SpeckleFilteringModule>("Speckle","SAR/Despeckle image");
//  model->RegisterModule<otb::SarIntensityModule>("SarIntensity","SAR/Compute intensity and log-intensity");
  
//  model->RegisterModule<otb::ChangeDetectionModule>("ChangeDetection","Filtering/Change Detection");
  
//  model->RegisterModule<otb::ProjectionModule>("Projection","Geometry/Reproject image");
//  model->RegisterModule<otb::SuperimpositionModule>("Superimposition","Geometry/Superimpose two images");
//  model->RegisterModule<otb::HomologousPointExtractionModule>("HomologousPoints", "Geometry/Homologous points extraction");
  
  // Module from the RT
//  model->RegisterModule<otb::VectorDataTransformModule>("VectorData Transform", "Geometry/VectorData  transform");
  

  // Launch Monteverdi
  view->InitWidgets();
  view->Show();

  Fl::lock();

  if (parseResult->IsOptionPresent("--InputList"))
    {
      int numberOfImage = parseResult->GetNumberOfParameters("--InputList");
      for (int i = 0; i < numberOfImage; i++)
        {
        Fl::check();
        std::vector<std::string> moduleVector;

        std::string filename = parseResult->GetParameterString("--InputList", i);
        if( itksys::SystemTools::FileExists(filename.c_str()) )
          {
            // Get the ModuleInstanceId
            std::string readerId = model->CreateModuleByKey("Reader");
            
            // Get the module itself
            otb::Module::Pointer module = model->GetModuleByInstanceId(readerId);
            
            // Simulate file chooser and ok callback
            otb::ReaderModule::Pointer readerModule = static_cast<otb::ReaderModule::Pointer>(dynamic_cast<otb::ReaderModule *>(module.GetPointer()));
            readerModule->vFilePath->value(parseResult->GetParameterString("--InputList", i).c_str());
            readerModule->Analyse();
            readerModule->bOk->do_callback();
            Fl::check();
          }
        else
          {
            itk::OStringStream oss;
            oss << "The file "<< filename <<" does not exist.";
            otb::MsgReporter::GetInstance()->SendError( oss.str().c_str() );
          }
        }
    }


  return Fl::run();
}
