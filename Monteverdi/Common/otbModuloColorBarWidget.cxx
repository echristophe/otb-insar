/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
See OTBCopyright.txt for details.


    This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE,  See the above copyright notices for more information.

=========================================================================*/
#include "otbModuloColorBarWidget.h"
#include "otbMath.h"
#include "FL/fl_draw.H"
#include <FL/Fl_Image.H>

namespace otb
{
/**
 * Constructor.
 */
ModuloColorBarWidget::ModuloColorBarWidget()
    : Fl_Window(0, 0, 0, 0, 0),m_OffsetValue(0.0)
{
  m_ScalarImage = SingleImageType::New();

}
/**
 * Destructor.
 */
ModuloColorBarWidget::~ModuloColorBarWidget()
{
}
/**
 * Show The widget.
 */
void ModuloColorBarWidget::Init(int x, int y, int w, int h, const char * l)
{

  this->label(l);
  this->resize(x, y, w, h);

  SingleImageType::RegionType region;
  SingleImageType::IndexType index;
  SingleImageType::SizeType size;

  size[0] = w;
  size[1] = h;
  region.SetSize(size);

  index[0] = 0;
  index[1] = 0;
  region.SetIndex(index);
  m_ScalarImage->SetRegions(region);
  m_ScalarImage->Allocate();

  typedef itk::ImageRegionIterator<SingleImageType> ImageRegionIteratorType;

  ImageRegionIteratorType it(m_ScalarImage, region);

  it.IsAtBegin();
  while (!it.IsAtEnd())
    {
    SingleImageType::IndexType index;
    index = it.GetIndex();
    SingleImageType::ValueType value;
    double tmp = static_cast<double> (255. * index[0] / (w-1) );
    value = static_cast<SingleImageType::ValueType> (tmp);
    it.Set(value);
    ++it;
    }
}


/**
 * Draw the widget
 */
void ModuloColorBarWidget::draw(void)
{

  ColorMapFilterType::Pointer  colorBarFilter = ColorMapFilterType::New();
  ModuloImageFilterType::Pointer moduloImage = ModuloImageFilterType::New();

  moduloImage->SetInput(m_ScalarImage);
  moduloImage->SetFactor(255.0);
  moduloImage->SetOffset(m_OffsetValue);

  colorBarFilter->SetInput(moduloImage->GetOutput());
  colorBarFilter->SetColormap(m_Colormap);
  colorBarFilter->Update();

  unsigned int height = colorBarFilter->GetOutput()->GetLargestPossibleRegion().GetSize(1);
  unsigned int width  = colorBarFilter->GetOutput()->GetLargestPossibleRegion().GetSize(0);

  unsigned int dim = width * height *3;

  typedef itk::ImageRegionConstIterator<ColorMapFilterType::OutputImageType> ImageRegionConstIteratorType;

  itk::VariableLengthVector<unsigned char>  colormap;
  colormap.SetSize(dim);

  SingleImageType::RegionType region;
  region = colorBarFilter->GetOutput()->GetRequestedRegion();

  ImageRegionConstIteratorType itColorBar(colorBarFilter->GetOutput(), region);

  unsigned int pos = 0;
  itColorBar.IsAtBegin();
  while (!itColorBar.IsAtEnd())
    {
      colormap[ pos    ] = itColorBar.Get()[0];
      colormap[ pos +1 ] = itColorBar.Get()[1];
      colormap[ pos +2 ] = itColorBar.Get()[2];
      pos +=3;
      ++itColorBar;
    }

  fl_draw_image(colormap.GetDataPointer(), 0 , 0, width, height, 3 );
}


}// end namespace otb
