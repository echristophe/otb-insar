/*=========================================================================

   Copyright 2012 Patrick IMBO
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
#ifndef __otbDynamicFringeViewerModule_cxx
#define __otbDynamicFringeViewerModule_cxx

#include <string>
#include <sstream>
#include <math.h>
#include "otbDynamicFringeViewerModule.h"
#include <FLU/Flu_File_Chooser.h>
#include <FL/Fl_Color_Chooser.H>
#include <FL/Fl_Text_Buffer.H>
#include "otbFltkFilterWatcher.h"
#include "otbMsgReporter.h"
#include "otbViewerConst.h"
#include <boost/algorithm/string.hpp>

namespace otb
{

/**
 * Constructor
 */
DynamicFringeViewerModule::DynamicFringeViewerModule() :
      m_ImageRenderingModel(), m_ImageView(), m_ImageController(), m_PixelDescriptionModel(), 
      m_PixelDescriptionView(), m_CurveWidget(), m_WindowsLayout(PACKED_WINDOWS_LAYOUT),
      m_DisplayMode(SLIDESHOW_DISPLAY_MODE),m_RenderingFunction(),
	  m_ColorMapRenderingModel(), m_ColorMapView(), m_ColorMapController(), m_FringeSpeed(0.0),
	  m_ColorMapOffset(0)
{
  // This module needs pipeline locking
  this->NeedsPipelineLockingOn();

  // Build rendering model
  m_ImageRenderingModel = RenderingModelType::New();
  m_ColorMapRenderingModel = RenderingModelType::New();

  m_PixelDescriptionModel = PixelDescriptionModelType::New();
  m_PixelDescriptionModel->SetLayers(m_ImageRenderingModel->GetLayers());

  // Build a view
  m_ImageView = ViewType::New();
  m_ColorMapView = ViewType::New();
  m_PixelDescriptionView = PixelDescriptionViewType::New();

  // Build a controller
  m_ImageController = WidgetControllerType::New();
  m_ColorMapController = WidgetControllerType::New();

  // Build the curve widget
  m_CurveWidget = CurvesWidgetType::New();
  m_CurveWidget->SetXAxisLabel(otbGetTextMacro("Pixels"));
  m_CurveWidget->SetYAxisLabel(otbGetTextMacro("Frequency"));

  // Wire the MVC
  m_ImageView->SetModel(m_ImageRenderingModel);
  m_ImageView->SetController(m_ImageController);
  m_PixelDescriptionView->SetModel(m_PixelDescriptionModel);

  m_ColorMapView->SetModel(m_ColorMapRenderingModel);
  m_ColorMapView->SetController(m_ColorMapController);

  // Add the resizing handler
  ResizingHandlerType::Pointer resizingHandler = ResizingHandlerType::New();
  resizingHandler->SetModel(m_ImageRenderingModel);
  resizingHandler->SetView(m_ImageView);
  m_ImageController->AddActionHandler(resizingHandler);

  // Add the change extract region handler
  ChangeRegionHandlerType::Pointer changeHandler = ChangeRegionHandlerType::New();
  changeHandler->SetModel(m_ImageRenderingModel);
  changeHandler->SetView(m_ImageView);
  m_ImageController->AddActionHandler(changeHandler);

  // Add the change scaled region handler
  ChangeScaledRegionHandlerType::Pointer changeScaledHandler = ChangeScaledRegionHandlerType::New();
  changeScaledHandler->SetModel(m_ImageRenderingModel);
  changeScaledHandler->SetView(m_ImageView);
  m_ImageController->AddActionHandler(changeScaledHandler);

  // Add the change scaled handler
  ChangeScaleHandlerType::Pointer changeScaleHandler = ChangeScaleHandlerType::New();
  changeScaleHandler->SetModel(m_ImageRenderingModel);
  changeScaleHandler->SetView(m_ImageView);
  m_ImageController->AddActionHandler(changeScaleHandler);

  // Add the pixel description action handler
  PixelDescriptionActionHandlerType::Pointer pixelActionHandler = PixelDescriptionActionHandlerType::New();
  pixelActionHandler->SetView(m_ImageView);
  pixelActionHandler->SetModel(m_PixelDescriptionModel);
  m_ImageController->AddActionHandler(pixelActionHandler);

  // Add the action handler for the arrow key
  ArrowKeyMoveActionHandlerType::Pointer arrowKeyMoveHandler = ArrowKeyMoveActionHandlerType::New();
  arrowKeyMoveHandler->SetModel(m_ImageRenderingModel);
  arrowKeyMoveHandler->SetView(m_ImageView);
  m_ImageController->AddActionHandler(arrowKeyMoveHandler);

  // Add the full widget handler
  DragFullActionHandlerType::Pointer dragHandler = DragFullActionHandlerType::New();
  dragHandler->SetModel(m_ImageRenderingModel);
  dragHandler->SetView(m_ImageView);
  m_ImageController->AddActionHandler(dragHandler);


  // Managed windows layout
  m_SplittedWindows = SplittedWidgetManagerType::New();
  m_PackedWindows = PackedWidgetManagerType::New();

  if (m_WindowsLayout == SPLITTED_WINDOWS_LAYOUT)
    {
    m_DisplayWindow = m_SplittedWindows;
    }
  else
    {
    m_DisplayWindow = m_PackedWindows;
    }
  m_DisplayWindow->RegisterFullWidget(m_ImageView->GetFullWidget());
  m_DisplayWindow->RegisterScrollWidget(m_ImageView->GetScrollWidget());
  m_DisplayWindow->RegisterZoomWidget(m_ImageView->GetZoomWidget());
  m_DisplayWindow->RegisterPixelDescriptionWidget(m_PixelDescriptionView->GetPixelDescriptionWidget());
  m_DisplayWindow->RegisterHistogramWidget(m_CurveWidget);

  // Data List Instance
  m_InputImage = SingleImageType::New();
  m_InputImageLayer = ImageLayerType::New();

  // Describe inputs
  this->AddInputDescriptor<SingleImageType> ("InputImage", otbGetTextMacro("Fringe image to display"), false, false);
  //this->AddTypeToInputDescriptor<SingleImageType> ("InputImage");
  //this->AddTypeToInputDescriptor<LabeledImageType> ("InputImage");
  //this->AddTypeToInputDescriptor<FloatImageWithQuicklook> ("InputImage");


  m_ColorMapFilter = ColorMapFilterType::New();
  m_ColorBarWidget = ColorBarWidgetType::New();
  m_ColorBarWidget->SetColormap(m_ColorMapFilter->GetColormap());

  m_ColorMapImageLayer = ImageLayerType::New();

  m_RenderingColorBarFilter = ColorMapFilterType::New();
  m_RenderingColorBarWidget = ColorBarWidgetType::New();
  m_RenderingColorBarWidget->SetColormap(m_RenderingColorBarFilter->GetColormap());

  // Build GUI
  this->Build();

  m_ColorBarWidget->Init(oColorBar->x(), oColorBar->y(),
                   oColorBar->w(), oColorBar->h(),"Color Bar");
  oColorBar->add(m_ColorBarWidget);
  oColorBar->box(FL_NO_BOX);
  m_ColorBarWidget->show();
  m_ColorBarWidget->redraw();

  m_RenderingColorBarWidget->Init(oRenderingColorBar->x(), oRenderingColorBar->y(),
                   oRenderingColorBar->w(), oRenderingColorBar->h(),"Color Bar");
  oRenderingColorBar->add(m_RenderingColorBarWidget);
  oRenderingColorBar->box(FL_NO_BOX);
  m_RenderingColorBarWidget->show();
  m_RenderingColorBarWidget->redraw();

  this->UpdateColorBar();

  // build the Screen shot GUI
  this->BuildScreenShot();
}

/**
 * Destructor
 */
DynamicFringeViewerModule::~DynamicFringeViewerModule()
{
}

/**
 * PrintSelf method
 */
void DynamicFringeViewerModule::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  // Call superclass implementation
  Superclass::PrintSelf(os, indent);
}

/**
 * The custom run command
 */
void DynamicFringeViewerModule::Run()
{
  // While the viewer is shown, it is busy
  this->BusyOn();


  m_InputImage = this->GetInputData<SingleImageType>("InputImage");
  if (m_InputImage.IsNull())
    {
      this->BusyOff();
      itkExceptionMacro(<< "Input image is NULL");
    }

  this->ColorMappingProcess();

  // Show everything
  m_DisplayWindow->Show();

  try
    {
    // Update the rendering model
    m_ImageRenderingModel->Update();
    m_ColorMapRenderingModel->Update();

    // Show the interface setup
    bSetupWindow->show();
    }
  catch (itk::ExceptionObject & err)
    {
    itk::OStringStream oss;
    oss << "Problem occured while loading input image. The following error was returned:\n";
    oss << err.GetDescription();
    MsgReporter::GetInstance()->SendError(oss.str());
    this->Quit();
    }
}



void 
DynamicFringeViewerModule
::UpdateColorBar()
{
  ColorMapFilterType::ColormapEnumType     colormapEnum;
  // By default colormapEnum is initialized at Red but following code
  // will used this variable only if its modified by the following switch
  // statement
  colormapEnum = ColorMapFilterType::Red;

  switch (iColorMap->value())
  {
    case 0 :
      m_ColormapName = "Red";
      colormapEnum = ColorMapFilterType::Red;
      break;
    case 1 :
      m_ColormapName = "Green";
      colormapEnum = ColorMapFilterType::Green;
      break;
    case 2 :
      m_ColormapName = "Blue";
      colormapEnum = ColorMapFilterType::Blue;
      break;
    case 3 :
      m_ColormapName = "Grey";
      colormapEnum = ColorMapFilterType::Grey;
      break;
    case 4 :
      m_ColormapName = "Hot";
      colormapEnum = ColorMapFilterType::Hot;
      break;
    case 5 :
      m_ColormapName = "Cool";
      colormapEnum = ColorMapFilterType::Cool;
      break;
    case 6 :
      m_ColormapName = "Spring";
      colormapEnum = ColorMapFilterType::Spring;
      break;
    case 7 :
      m_ColormapName = "Summer";
      colormapEnum = ColorMapFilterType::Summer;
      break;
    case 8 :
      m_ColormapName = "Autumn";
      colormapEnum = ColorMapFilterType::Autumn;
      break;
    case 9 :
      m_ColormapName = "Winter";
      colormapEnum = ColorMapFilterType::Winter;
      break;
    case 10 :
      m_ColormapName = "Copper";
      colormapEnum = ColorMapFilterType::Copper;
      break;
    case 11 :
      m_ColormapName = "Jet";
      colormapEnum = ColorMapFilterType::Jet;
      break;
    case 12 :
      m_ColormapName = "HSV";
      colormapEnum = ColorMapFilterType::HSV;
      break;
    case 13 :
      m_ColormapName = "OverUnder";
      colormapEnum = ColorMapFilterType::OverUnder;
      break;
      break;
    default:
      itkExceptionMacro(<< "Colormap not implemented");
      break;
  }

  m_ColorMapFilter->SetColormap(colormapEnum);
  m_ColorBarWidget->SetColormap(m_ColorMapFilter->GetColormap());
  m_ColorBarWidget->redraw();

  m_RenderingColorBarFilter->SetColormap(colormapEnum);
  m_RenderingColorBarWidget->SetColormap(m_RenderingColorBarFilter->GetColormap());
  m_RenderingColorBarWidget->redraw();

}

/**
 *
 */
void DynamicFringeViewerModule::ColorMappingProcess()
{

	m_RenderingColorBarWidget->SetOffsetValue(m_ColorMapOffset);
	m_RenderingColorBarWidget->redraw();

  
	ModuloImageFilterType::Pointer moduloImage = ModuloImageFilterType::New();

	moduloImage->SetInput(m_InputImage);
	moduloImage->SetFactor(255.0);
	moduloImage->SetOffset(m_ColorMapOffset);

	ColorMapFilterType::Pointer colorMapImage = ColorMapFilterType::New();
	colorMapImage->SetInput(moduloImage->GetOutput());
    colorMapImage->SetColormap(m_ColorMapFilter->GetColormap());

	RGBtoVectorImageCastFilterType::Pointer RGBtoVectorImage = RGBtoVectorImageCastFilterType::New();
	RGBtoVectorImage->SetInput(colorMapImage->GetOutput());

	// Generate the layer
	ImageLayerGeneratorPointerType generator = ImageLayerGeneratorType::New();
	generator->SetImage(RGBtoVectorImage->GetOutput());

	try
		  {
		FltkFilterWatcher qlwatcher(generator->GetProgressSource(), 0, 0, 200, 20,
			                          otbGetTextMacro("Generating QuickLook ..."));
		generator->GenerateLayer();
      }
    catch (itk::ExceptionObject & err)
      {
      itk::OStringStream oss;
      oss << "Problem occurred while generation of QuickLook. The following error was returned:\n";
      oss << err.GetDescription();
      MsgReporter::GetInstance()->SendError(oss.str());
      this->Quit();
      }
    m_InputImageLayer = generator->GetLayer();

    m_RenderingFunction = generator->GetRenderingFunction();


    // Add the generated layer to the rendering model
    m_ImageRenderingModel->AddLayer(generator->GetLayer());

    // Show everything
    m_DisplayWindow->Show();

    // Update the rendering model
    m_ImageRenderingModel->Update();

    m_InputImageLayer->SetVisible(true);
}



/**
 *
 */
void DynamicFringeViewerModule::ActivateSlideShowMode()
{
  if (m_DisplayMode == SLIDESHOW_DISPLAY_MODE)
    {
    return;
    }

  m_DisplayMode = SLIDESHOW_DISPLAY_MODE;

  // Activate slide show related widgets
  bPreviousImage->activate();
  bNextImage->activate();

  ShowSelectedImages();
}


/**
 *
 */
void DynamicFringeViewerModule::ShowPreviousImage()
{
	m_ColorMapOffset -= m_FringeSpeed ;
	m_ColorMapOffset =  fmod(m_ColorMapOffset,255.0) ;
	ShowSelectedImages();
}

/**
 *
 */
void DynamicFringeViewerModule::ShowNextImage()
{
	m_ColorMapOffset += +m_FringeSpeed ;
	m_ColorMapOffset =  fmod(m_ColorMapOffset,255.0) ;

	ShowSelectedImages();
}

/**
 *
 */
void DynamicFringeViewerModule::ShowSelectedImages()
{
  ImageLayerType::Pointer imageLayer;
  m_ImageRenderingModel->ClearLayers();

  this->ColorMappingProcess();


  m_ColorMapRenderingModel->AddLayer(m_ColorMapImageLayer);
  
  const DataObjectWrapper& dow = this->GetInputDataDescriptorByKey(std::string("InputImage")).GetNthData(0);
  std::ostringstream title;
  title << "[" << dow.GetSourceInstanceId() << "] " << dow.GetSourceOutputKey();
  bSetupWindow->copy_label(title.str().c_str());
  m_SplittedWindows->SetLabel(title.str().c_str());
  m_PackedWindows->SetLabel(title.str().c_str());

  m_ImageRenderingModel->Update();
  m_ColorMapRenderingModel->Update();
  
  RedrawWidget();
}






/**
 *
 */
void DynamicFringeViewerModule::RedrawWidget()
{
  m_DisplayWindow->Refresh();
}




/**
 * Method call when click Apply button in contrast stretch gui.
 */
void DynamicFringeViewerModule::ApplyModuloRendering()
{
	  m_ImageRenderingModel->Update();
	  m_ColorMapRenderingModel->Update();

}



/**
 * Display the viewed elements in splitted windows
 */
void DynamicFringeViewerModule::SplittedLayout()
{
  if (m_WindowsLayout == SPLITTED_WINDOWS_LAYOUT)
    {
    return;
    }

  UpdateWindowsLayout(SPLITTED_WINDOWS_LAYOUT);
}

/**
 * Display the viewed elements in packed windows
 */
void DynamicFringeViewerModule::PackedLayout()
{
  if (m_WindowsLayout == PACKED_WINDOWS_LAYOUT)
    {
    return;
    }

  UpdateWindowsLayout(PACKED_WINDOWS_LAYOUT);
}

/**
 * Set the display mode to "windowsLayout" value
 */
void DynamicFringeViewerModule::UpdateWindowsLayout(const WindowsLayoutEnumType windowsLayout)
{
  WidgetManagerPointerType tmpHandler = m_DisplayWindow;

  m_DisplayWindow->Hide();
  m_DisplayWindow->UnRegisterAll();

  if (windowsLayout == SPLITTED_WINDOWS_LAYOUT)
    {
    m_DisplayWindow = m_SplittedWindows;
    }
  else
    {
    m_DisplayWindow = m_PackedWindows;
    }

  m_DisplayWindow->RegisterFullWidget(m_ImageView->GetFullWidget());
  m_DisplayWindow->RegisterScrollWidget(m_ImageView->GetScrollWidget());
  m_DisplayWindow->RegisterZoomWidget(m_ImageView->GetZoomWidget());
  m_DisplayWindow->RegisterPixelDescriptionWidget(m_PixelDescriptionView->GetPixelDescriptionWidget());
  m_DisplayWindow->RegisterHistogramWidget(m_CurveWidget);

  tmpHandler->Hide();
  tmpHandler->UnRegisterAll();

  m_DisplayWindow->Refresh();
  m_DisplayWindow->Show();

  m_WindowsLayout = windowsLayout;
}


/**
 *
 */
void DynamicFringeViewerModule::ShowHide()
{
  m_DisplayWindow->Show();
}

/**
 *
 */
void DynamicFringeViewerModule::Quit()
{
  this->Hide();
  // Once module is closed, it is no longer busy
  this->BusyOff();
}


/**
 * Proceed to screen shot
 */
void DynamicFringeViewerModule::ScreenShot()
{
  if (rbScreenZoom->value() == 1)
    {
    this->SaveScreenShot("Save Zoom view screen shot as...", m_ImageView->GetZoomWidget());
    }

  if (rbScreenFull->value() == 1)
    {
    this->SaveScreenShot("Save Full view screen shot as...", m_ImageView->GetFullWidget());
    }

  if (rbScreenNav->value() == 1)
    {
    this->SaveScreenShot("Save Navigation view screen shot as...", m_ImageView->GetScrollWidget());
    }

  if (rbScreenNav->value() != 1 && rbScreenFull->value() != 1 && rbScreenZoom->value() != 1)
    {
    MsgReporter::GetInstance()->SendError("No view selected for screen shot...");
    }

  wScreenShot->hide();
}


void DynamicFringeViewerModule::SaveScreenShot(const char * winLab, ImageWidgetType * widget)
{
  const char * filename = NULL;

  ScreenShotFilterType::Pointer screener = ScreenShotFilterType::New();
  screener->SetNumberOfChannels(3);
  screener->SetInverseXSpacing(true);

  filename = flu_file_chooser(otbGetTextMacro(winLab), "*.*", "");

  if (filename == NULL)
    {
    MsgReporter::GetInstance()->SendError("Empty file name!");
    return;
    }

  screener->SetFileName(filename);
  screener->SetBuffer(widget->GetOpenGlBuffer());
  screener->SetImageSize(widget->GetOpenGlBufferedRegion().GetSize());
  screener->Update();
}


/**
 *
 */
void DynamicFringeViewerModule::Hide()
{
  // Hide the main window
  m_DisplayWindow->Hide();
  // Hide the Screen shot Window
  wScreenShot->hide();
  // Hide the Setup Propreties Window
  bSetupWindow->hide();
}

/**
 * Proceed to screen shot
 */
void DynamicFringeViewerModule::SetFringeSpeed()
{
	m_FringeSpeed = guiFringeSpeed->value();
	this->RedrawWidget();
}


} // End namespace otb

#endif
