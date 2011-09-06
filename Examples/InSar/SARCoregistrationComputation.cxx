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
#include "otbStreamingResampleImageFilter.h"

#include "itkImageRegistrationMethod.h"
#include "otbGenericRSTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "otbComplexInterpolateImageFunction.h"
#include "otbWindowedSincInterpolateImageBlackmanFunction.h"

#include "itkGradientDescentOptimizer.h"
#include "otbLeeImageFilter.h"

#include "itkResampleImageFilter.h"

#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "otbComplexToIntensityImageFilter.h"

// Command line:

// ./SARCoregistrationComputation ...


/** Quick function for rounding */
int roundValue(double value)
{
  return ((value < 0.0) ?
      -std::ceil( std::abs(value) -0.5 )
     : std::ceil( std::abs(value) -0.5 )
         );
}

int main(int argc, char* argv[])
{

  if (argc != 4)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile";
    std::cerr << " registeredSlaveImageFile" << std::endl;
    return EXIT_FAILURE;
    }

  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,2> ImageType;

  typedef otb::ImageFileReader<ImageType> ReaderType;
  typedef otb::StreamingImageFileWriter<ImageType> WriterType;

  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(argv[1]);
  slave->SetFileName(argv[2]);

  master->UpdateOutputInformation();
  slave->UpdateOutputInformation();

  std::cout << "Master size: " << std::endl;
  std::cout << master->GetOutput()->GetLargestPossibleRegion() << std::endl;

  std::cout << "Slave size: " << std::endl;
  std::cout << slave->GetOutput()->GetLargestPossibleRegion() << std::endl;
  
  /* Set up registration framework */
  typedef otb::GenericRSTransform< double, 2 >										TransformType;
  typedef itk::GradientDescentOptimizer												OptimizerType;
  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType >		MetricType;

  typedef otb::Function::BlackmanWindowFunction<double>			FunctionType;
  typedef itk::ConstantBoundaryCondition<ImageType>				BoundaryConditionType;
  typedef double												CoordRepType;

  typedef otb::ComplexInterpolateImageFunction< ImageType, FunctionType, BoundaryConditionType, CoordRepType >			InterpolatorType;
  typedef itk::ImageRegistrationMethod< ImageType, ImageType >						RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform(transform);
  registration->SetInterpolator(interpolator);

  registration->SetFixedImage(master->GetOutput());
  registration->SetMovingImage(slave->GetOutput());


  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y

  registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetLearningRate( 70.0 );
  optimizer->SetNumberOfIterations( 100 );

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  //scale values for amplitude
//   optimizerScales[0]= 10000000000;
//   optimizerScales[1]= 10000000000;
//  optimizerScales[0]= 100000;// related to the amplitude ERS
//  optimizerScales[1]= 20000;
  optimizerScales[0]= 100000*1000000L;// related to the amplitude Palsar
  optimizerScales[1]= 20000*1000000L;
  optimizer->SetScales(optimizerScales);

  registration->Update();

  /*
  //second part of the registration
  registration->SetInitialTransformParameters( registration->GetLastTransformParameters());
  optimizer->SetLearningRate( 50.0 );
  optimizer->SetNumberOfIterations( 150 );
  registration->Update();
*/
  ParametersType finalParameters = registration->GetLastTransformParameters();

  std::cout << "The final estimated translation is: " << finalParameters << std::endl;

  finalParameters[0] = roundValue(finalParameters[0]);
  finalParameters[1] = roundValue(finalParameters[1]);

  std::cout << "The integer translation is: " << finalParameters << std::endl;

  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetParameters( finalParameters );

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( finalTransform );
  resample->SetInput( slave->GetOutput() );

  resample->SetSize(    master->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  master->GetOutput()->GetOrigin() );
  resample->SetOutputSpacing( master->GetOutput()->GetSpacing() );
  resample->SetDefaultPixelValue( 0 );

  typedef otb::ImageFileWriter<ImageType> WriterFixedType;
  WriterFixedType::Pointer writer = WriterFixedType::New();
  writer->SetFileName( "slave-registered.tif" );
  writer->SetInput( resample->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
