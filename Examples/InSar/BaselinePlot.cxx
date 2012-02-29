
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

// Command line:
//
// ./bin/BaselinePlot ~/project/Images/TSX1_SAR__SSC______HS_S_SRA_20090212T204239_20090212T204240/TSX1_SAR__SSC______HS_S_SRA_20090212T204239_20090212T204240.xml ~/project/Images/TSX1_SAR__SSC______HS_S_SRA_20090223T204240_20090223T204241/TSX1_SAR__SSC______HS_S_SRA_20090223T204240_20090223T204241.xml

#include <iomanip>

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbBaselineCalculator.h"
#include "otbLengthOrientationBaselineFunctor.h"
#include "otbPlatformPositionToBaselineCalculator.h"

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkAxis.h>
#include <vtkPlotPoints.h>

int main(int argc, char* argv[])
{

  if (argc != 3)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 2 ;

  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,Dimension> ImageType;

  typedef otb::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(argv[1]);
  slave->SetFileName(argv[2]);

  master->UpdateOutputInformation();
  slave->UpdateOutputInformation();
  typedef otb::Functor::LengthOrientationBaselineFunctor	BaselineFunctorType;
  typedef otb::BaselineCalculator<BaselineFunctorType>    BaselineCalculatorType;
  typedef BaselineCalculatorType::PlateformPositionToBaselineCalculatorType PlateformPositionToBaselineCalculatorType;

  BaselineCalculatorType::Pointer baselineCalculator = BaselineCalculatorType::New();

  BaselineCalculatorType::PlateformPositionToBaselinePointer plateformPositionToBaseline =  PlateformPositionToBaselineCalculatorType::New();
  plateformPositionToBaseline->SetMasterPlateform(master->GetOutput()->GetImageKeywordlist());
  plateformPositionToBaseline->SetSlavePlateform(slave->GetOutput()->GetImageKeywordlist());

  baselineCalculator->SetPlateformPositionToBaselineCalculator(plateformPositionToBaseline);
  baselineCalculator->Compute(otb::Functor::LengthOrientationBaselineFunctor::Length);

  std::vector<BaselineCalculatorType::PointType> pointImage;
  pointImage.clear();
  std::vector<double> baselineImage;
  baselineImage.clear();
  baselineCalculator->ExtractBaseline(otb::Functor::LengthOrientationBaselineFunctor::Length, 
						pointImage, baselineImage);


  /** Table extract from function*/
  // Create a table with some points in it
  vtkSmartPointer<vtkTable> tableFunction = vtkSmartPointer<vtkTable>::New();
 
  vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
  arrX->SetName("Line index");
  tableFunction->AddColumn(arrX);
 
  vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
  arrC->SetName("Baseline");
  tableFunction->AddColumn(arrC);
 
  // Fill in the table with some example values
  int numPoints = 100;

  float inc = master->GetOutput()->GetLargestPossibleRegion().GetSize(0) / (numPoints-1);
  tableFunction->SetNumberOfRows(numPoints);
  for (int i = 0; i < numPoints; ++i)
  {
    tableFunction->SetValue(i, 0, i * inc);
	tableFunction->SetValue(i, 1, baselineCalculator->EvaluateBaseline(i * inc,0.0));
  }
 
  /** Table extract from function*/
  // Create a table with some points in it
  vtkSmartPointer<vtkTable> tableScatter = vtkSmartPointer<vtkTable>::New();
 
  vtkSmartPointer<vtkFloatArray> scatterX = vtkSmartPointer<vtkFloatArray>::New();
  scatterX->SetName("Line index");
  tableScatter->AddColumn(scatterX);
 
  vtkSmartPointer<vtkFloatArray> scatterY = vtkSmartPointer<vtkFloatArray>::New();
  scatterY->SetName("Baseline");
  tableScatter->AddColumn(scatterY);
 
  // Fill in the table with some example values
  numPoints = baselineImage.size();
  tableScatter->SetNumberOfRows(numPoints);
  for (int i = 0; i < numPoints; ++i)
  {
    tableScatter->SetValue(i, 0, pointImage[i][0]);
	tableScatter->SetValue(i, 1, baselineImage[i]);
  }

  /***********************/

  // Set up the view
  vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
  view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
 
  // Add multiple line plots, setting the colors etc
  vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
  chart->GetAxis(0)->SetTitle("Baseline");
  chart->GetAxis(1)->SetTitle("Line index");

  view->GetScene()->AddItem(chart);
 
  vtkPlot *line = chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION <= 5
  line->SetInput(tableFunction, 0, 1);
#else
  line->SetInputData(tableFunction, 0, 1);
#endif
  line->SetColor(255, 0, 0, 255);
  line->SetWidth(1.0);
  line = chart->AddPlot(vtkChart::LINE);

  vtkPlot * points = chart->AddPlot(vtkChart::POINTS);
#if VTK_MAJOR_VERSION <= 5
  points->SetInput(tableScatter, 0, 1);
#else
  points->SetInputData(tableScatter, 0, 1);
#endif
  points->SetColor(0, 0, 0, 255);
  points->SetWidth(1.0);
  vtkPlotPoints::SafeDownCast(points)->SetMarkerStyle(vtkPlotPoints::PLUS);
 
  view->GetRenderWindow()->SetMultiSamples(0);
 
  // Start interactor
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();
 

  std::cout << "----------------------" << std::endl;
 
  return EXIT_SUCCESS;
}
