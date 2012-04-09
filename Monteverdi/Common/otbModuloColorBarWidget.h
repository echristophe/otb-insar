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
#ifndef __otbModuloColorBarWidget_h
#define __otbModuloColorBarWidget_h

#include "FL/Fl_Window.H"
#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkMacro.h"

#include "otbModule.h"
#include "otbImage.h"
#include "otbVectorImage.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "otbRGBImageToVectorImageCastFilter.h"
#include "otbModuloImageFilter.h"

namespace otb
{
/** \class ModuloColorBarWidget
 *  \brief
 *
 */
class ModuloColorBarWidget
      : public itk::Object, public Fl_Window
{
public:
  /** Standard class typedefs */
  typedef ModuloColorBarWidget                Self;
  typedef itk::Object                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Runtime information */
  itkTypeMacro(ModuloColorBarWidget, Object);

  // Convenient typedefs
  typedef TypeManager::Floating_Point_Precision     Floating_Point_PrecisionType;
  typedef TypeManager::Floating_Point_Image         SingleImageType;
  typedef TypeManager::Floating_Point_VectorImage   OutputImageType;


  typedef itk::RGBPixel<unsigned char>              RGBPixelType;
  typedef otb::Image<RGBPixelType, 2>               RGBImageType;

  // ColorMapping Class typedefs
  typedef otb::ModuloImageFilter<SingleImageType, SingleImageType>  ModuloImageFilterType;
  typedef itk::ScalarToRGBColormapImageFilter<SingleImageType, RGBImageType>  ColorMapFilterType;
  typedef ColorMapFilterType::ColormapType ColormapType;
  typedef otb::RGBImageToVectorImageCastFilter<RGBImageType, OutputImageType> RGBtoVectorImageCastFilterType;


  /** Initialize the widget */
  virtual void Init(int x, int y, int w, int h, const char * l);


  /** Set Colormap */
  itkSetMacro(Colormap, ColormapType::Pointer);

 itkSetClampMacro(OffsetValue,double,0.0,255.0);
 itkGetConstMacro(OffsetValue,double);

protected:
  /** Constructor */
  ModuloColorBarWidget();
  /** Destructor */
  ~ModuloColorBarWidget();

  /** Draw the widget */
  virtual void draw(void);

private:
  ModuloColorBarWidget(const Self&); // purposely not implemented
  void operator=(const Self&); // purposely not implemented

  ColormapType::Pointer     m_Colormap;
  SingleImageType::Pointer  m_ScalarImage;
  double                    m_OffsetValue; 
};
} // end namespace otb
#endif
