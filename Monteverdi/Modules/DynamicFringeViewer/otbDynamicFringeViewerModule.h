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
#ifndef __otbDynamicFringeViewerModule_h
#define __otbDynamicFringeViewerModule_h

// include the base class
#include "otbModule.h"

// include the OTB elements
#include "otbVectorImage.h"
#include "otbVectorData.h"

//standardImage Viewer
#include "otbDynamicFringeViewerModuleGUI.h"

#include "otbImage.h"
#include "otbVectorImage.h"
#include "itkRGBAPixel.h"
#include "otbImageLayer.h"
#include "otbImageLayerRenderingModel.h"
#include "otbImageView.h"
#include "otbImageWidgetController.h"
#include "otbImageLayerGenerator.h"
#include "otbWidgetResizingActionHandler.h"
#include "otbArrowKeyMoveActionHandler.h"
#include "otbChangeScaledExtractRegionActionHandler.h"
#include "otbChangeExtractRegionActionHandler.h"
#include "otbChangeScaleActionHandler.h"
#include "otbHistogramActionHandler.h"
#include "otbCurves2DWidget.h"
#include "otbHistogramCurve.h"
#include "otbVerticalAsymptoteCurve.h"
#include "otbPixelDescriptionModel.h"
#include "otbPixelDescriptionActionHandler.h"
#include "otbPixelDescriptionView.h"
#include "otbStandardRenderingFunction.h"
#include "otbNoStretchRenderingFunction.h"
#include "otbSquareRootRenderingFunction.h"
#include "otbGaussianRenderingFunction.h"
#include "otbImageToVectorImageCastFilter.h"
#include "otbDragFullWindowActionHandler.h"
#include "otbScalarBufferToImageFileWriter.h"

#include "otbAmplitudeFunctor.h"
#include "otbPhaseFunctor.h"

#include "otbWidgetManager.h"
#include "otbPackedWidgetManager.h"
#include "otbSplittedWidgetManager.h"
#include "otbMonteverdiEnum.h"
#include "otbUniformAlphaBlendingFunction.h"
#include "otbBlendingFunction.h"
#include "itkRescaleIntensityImageFilter.h"

#include "otbModuloColorBarWidget.h"
#include "otbModuloImageFilter.h"
#include "otbRGBImageToVectorImageCastFilter.h"
#include "otbVectorImage.h"


namespace otb
{


/** \class DynamicFringeViewerModule
 *  \brief
 *
 *  \sa DataObjectWrapper, DataDescriptor, DataDescriptor
 */
class ITK_EXPORT DynamicFringeViewerModule
  : public Module, public DynamicFringeViewerModuleGUI
{
public:
  /** Standard class typedefs */
  typedef DynamicFringeViewerModule                  Self;
  typedef Module                        Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** New macro */
  itkNewMacro(Self);

  /** Type macro */
  itkTypeMacro(DynamicFringeViewerModule, Module);

  /** Dataset */
  typedef TypeManager::Floating_Point_Precision     PixelType;
  typedef TypeManager::Floating_Point_Image         SingleImageType;
  typedef TypeManager::Floating_Point_VectorImage   VectorImageType;

  /** WidgetManager related types definition */
  typedef WidgetManager                             WidgetManagerType;
  typedef WidgetManagerType::Pointer                WidgetManagerPointerType;

  typedef PackedWidgetManager                       PackedWidgetManagerType;
  typedef PackedWidgetManagerType::Pointer          PackedWidgetManagerPointerType;

  typedef SplittedWidgetManager                     SplittedWidgetManagerType;
  typedef SplittedWidgetManagerType::Pointer        SplittedWidgetManagerPointerType;

  /** Input image convenient typedefs */
  typedef VectorImageType::PixelType                      InputPixelType;
  typedef VectorImageType::InternalPixelType              InputInternalPixelType;

  /** Output image type */
  typedef itk::RGBPixel<unsigned char>              RGBPixelType;
  typedef Image<RGBPixelType, 2>                    OutputImageType;

  /** Image layer type */
  typedef ImageLayer<VectorImageType, OutputImageType>    ImageLayerType;
  typedef ImageLayerType::Pointer                   ImageLayerPointerType;
  typedef ImageLayerType::HistogramType             HistogramType;

  /** Image layer generator type */
  typedef ImageLayerGenerator<ImageLayerType>       ImageLayerGeneratorType;
  typedef ImageLayerGeneratorType::Pointer          ImageLayerGeneratorPointerType;

  /** Rendering model type */
  typedef ImageLayerRenderingModel<OutputImageType> RenderingModelType;
  typedef RenderingModelType::Pointer               RenderingModelPointerType;

  /** View type */
  typedef ImageView<RenderingModelType>             ViewType;
  typedef ViewType::Pointer                         ViewPointerType;
  typedef ViewType::ImageWidgetType                 ImageWidgetType;

  /** Widget controller */
  typedef ImageWidgetController                     WidgetControllerType;
  typedef WidgetControllerType::Pointer             WidgetControllerPointerType;

  /** Curves 2D widget */
  typedef Curves2DWidget                            CurvesWidgetType;
  typedef CurvesWidgetType::Pointer                 CurvesWidgetPointerType;

  /** Rendering function */
  typedef ImageLayerGeneratorType::RenderingFunctionType StandardRenderingFunctionType;
  typedef StandardRenderingFunctionType::Pointer         StandardRenderingFunctionPointerType;


  /** Rendering function */
  typedef Function::RenderingFunction<SingleImageType::PixelType, RGBPixelType> RenderingFunctionType;

  typedef otb::Function::ChannelSelectorFunctor<PixelType> PixelRepresentationFunctionType;
  typedef otb::Function::Modulo<PixelType, PixelType>       ModuloFunctionType;

  typedef Function::StandardRenderingFunction<SingleImageType::PixelType,
                                              RGBPixelType,
											  PixelRepresentationFunctionType,
											  ModuloFunctionType				
											  > ModuloRenderingFunctionType;
  typedef ModuloRenderingFunctionType::Pointer  ModuloRenderingFunctionPointerType;
  /** Standard action handlers */
  typedef otb::WidgetResizingActionHandler
  <RenderingModelType, ViewType>                     ResizingHandlerType;
  typedef otb::ChangeScaledExtractRegionActionHandler
  <RenderingModelType, ViewType>                     ChangeScaledRegionHandlerType;
  typedef otb::ChangeExtractRegionActionHandler
  <RenderingModelType, ViewType>                     ChangeRegionHandlerType;
  typedef otb::ChangeScaleActionHandler
  <RenderingModelType, ViewType>                     ChangeScaleHandlerType;
  typedef otb::ArrowKeyMoveActionHandler
  <RenderingModelType, ViewType>                     ArrowKeyMoveActionHandlerType;
  typedef otb::HistogramActionHandler
    <RenderingModelType, ViewType,
     RenderingFunctionType>                         HistogramActionHandlerType;
  typedef otb::DragFullWindowActionHandler
  <RenderingModelType, ViewType>                    DragFullActionHandlerType;

  /** Pixel description */
  typedef PixelDescriptionModel<OutputImageType>     PixelDescriptionModelType;
  typedef  PixelDescriptionModelType::Pointer        PixelDescriptionModelPointerType;
  typedef PixelDescriptionActionHandler
    < PixelDescriptionModelType, ViewType>           PixelDescriptionActionHandlerType;
  typedef otb::PixelDescriptionView
    < PixelDescriptionModelType >                    PixelDescriptionViewType;
  typedef  PixelDescriptionViewType::Pointer         PixelDescriptionViewPointerType;

  /** Cast SingleImage to VectorImageType*/
  typedef ImageToVectorImageCastFilter<SingleImageType, VectorImageType>     CastSingleImageFilter;

  typedef otb::ModuloImageFilter<SingleImageType, SingleImageType>  ModuloImageFilterType;

  /** Screen shot typedefs */
  typedef ScalarBufferToImageFileWriter<unsigned char, unsigned char> ScreenShotFilterType;

  typedef itk::ScalarToRGBColormapImageFilter<SingleImageType, OutputImageType>              ColorMapFilterType;
  typedef otb::RGBImageToVectorImageCastFilter<OutputImageType, VectorImageType>             RGBtoVectorImageCastFilterType;
  typedef itk::ScalarToRGBColormapImageFilter<SingleImageType, OutputImageType>              ColorMapFilterType;
  //typedef otb::ColorBarWidget ColorBarWidgetType;
  typedef otb::ModuloColorBarWidget ColorBarWidgetType;

  /** Structure containing parameters of viewer setup */
  struct _ViewerSetupStruct
  {
    unsigned int pGrayChannel;
    unsigned int pRedChannel;
    unsigned int pGreenChannel;
    unsigned int pBlueChannel;
    VectorImageType::InternalPixelType pNoData;
    double pLowerQuantile;
    double pUpperQuantile;
    double pStandardDeviation;
    StretchResolutionEnumType pStretchResolution;
    bool pRGBMode;
  };

  typedef struct _ViewerSetupStruct ViewerSetupStructType;


 protected:
  /** Constructor */
  DynamicFringeViewerModule();

  /** Destructor */
  virtual ~DynamicFringeViewerModule();

  /** PrintSelf method */
  virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

  /** The custom run command */
  virtual void Run();

  /** CallBacks Implementations*/
  virtual void ActivateSlideShowMode();
  virtual void ShowPreviousImage();
  virtual void ShowNextImage();
  virtual void ShowSelectedImages();
  virtual void RedrawWidget();
  virtual void ScreenShot();

  /** Setup viewer Callbacks*/
  virtual void UpdateWindowsLayout(const WindowsLayoutEnumType windowsLayout);
  virtual void PackedLayout();
  virtual void SplittedLayout();
  virtual void ShowHide();
  virtual void Quit();
  virtual void SetFringeSpeed();
  virtual void ApplyModuloRendering();

  /** Hide every window */
  virtual void Hide();

  virtual void ColorMappingProcess();

  virtual void UpdateColorBar();
  void SaveScreenShot( const char * winLab, ImageWidgetType * widget);

private:
  DynamicFringeViewerModule(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  SingleImageType::Pointer                   m_InputImage;

  /** The image layers */
  ImageLayerType::Pointer                  m_InputImageLayer;
  ViewerSetupStructType                    m_ViewerSetupStruct;

   /** The Image rendering model */
  RenderingModelPointerType                m_ImageRenderingModel;

  /** The pixel description model */
  PixelDescriptionModelPointerType         m_PixelDescriptionModel;

  /** The view */
  ViewPointerType                          m_ImageView;

    /** The pixel description view */
  PixelDescriptionViewPointerType          m_PixelDescriptionView;

  /** Curve widget */
  CurvesWidgetPointerType                  m_CurveWidget;

  /** The widget controller */
  WidgetControllerPointerType              m_ImageController;

  /** Rendering function */
  StandardRenderingFunctionPointerType           m_RenderingFunction;

  /** Management of the windows layout (splitted or packed windows) */
  WidgetManagerPointerType                 m_DisplayWindow;
  PackedWidgetManagerPointerType           m_PackedWindows;
  SplittedWidgetManagerPointerType         m_SplittedWindows;
  WindowsLayoutEnumType                    m_WindowsLayout;

  /** Management of the display mode (image overlay with transparency or slide show) */
  DisplayModeEnumType                      m_DisplayMode;

  ColorBarWidgetType::Pointer              m_ColorBarWidget;
  ColorMapFilterType::Pointer              m_ColorMapFilter;
  std::string                              m_ColormapName;

  ColorBarWidgetType::Pointer              m_RenderingColorBarWidget;
  ColorMapFilterType::Pointer              m_RenderingColorBarFilter;
  double								   m_FringeSpeed;	

  /** ColorMap */
  ImageLayerType::Pointer                  m_ColorMapImageLayer;
  RenderingModelPointerType                m_ColorMapRenderingModel;
  ViewPointerType                          m_ColorMapView;
  WidgetControllerPointerType              m_ColorMapController;
  double								   m_ColorMapOffset;	

};


} // End namespace otb

#endif /* __otbDynamicFringeViewerModule_h */
