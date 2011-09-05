/*=========================================================================

   Copyright 2009 Emmanuel Christophe
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
#include "otbForwardSensorModel.h"
#include "otbInverseSensorModel.h"
#include "otbCompositeTransform.h"
#include "otbStreamingResampleImageFilter.h"
#include "otbExtractROI.h"
#include "itkUnaryFunctorImageFilter.h"

#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

// #include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkMeanImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "otbLeeImageFilter.h"

#include "otbBinaryFunctorNeighborhoodImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "otbAmplitudePhaseToRGBFunctor.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "otbComplexToIntensityImageFilter.h"

// Command line:

// ./BasicInterferogramComputation ~/data/ERS-Gamma/05721/05721.slc.hdr ~/data/ERS-Gamma/16242/16242.slc.hdr interferogram.hdr colorInterferogram.png

// ./BasicInterferogramComputation /home/christop/data2/Palsar/486_0010_20080510_FBD_11/VOL-ALPSRP122180010-H1.1__A /home/christop/data2/Palsar/486_0010_20090928_FBD_11/VOL-ALPSRP195990010-H1.1__A interferogramPalsar.hdr colorInterferogramPalsar.png

/** Helper class to display progress of the registration */
class CommandIterationUpdate : public itk::Command
{
  public:
    typedef  CommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer;
    itkNewMacro( Self );

  protected:
    CommandIterationUpdate() {};

  public:

//     typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
    typedef itk::GradientDescentOptimizer     OptimizerType;
    typedef const OptimizerType                         *OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
          dynamic_cast< OptimizerPointer >( object );

      if ( ! itk::IterationEvent().CheckEvent( &event ) )
      {
        return;
      }

      std::cout << optimizer->GetCurrentIteration() << " = ";
      std::cout << optimizer->GetValue() << " : ";
      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }

};

/**
 * Helper class to record progress of the registration
 * the curves can be visualized with the plotRegistration.py
*/
class CommandIterationUpdatePlot : public itk::Command
{
  public:
    typedef  CommandIterationUpdatePlot   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer;
    itkNewMacro( Self );

    void SetFileName(std::string filename)
    {
      m_filename=filename;
      //remove file
      std::ofstream file;
      file.open(m_filename.c_str(),std::ios::trunc);
      file.close();
    }

  protected:
    CommandIterationUpdatePlot()
    {
      m_filename = "out-GitFobquen3.dat";
    };

    std::string m_filename;


  public:

//     typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
    typedef itk::GradientDescentOptimizer     OptimizerType;
    typedef const OptimizerType                         *OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
          dynamic_cast< OptimizerPointer >( object );

      if ( ! itk::IterationEvent().CheckEvent( &event ) )
      {
        return;
      }
      std::ofstream file;
      file.open(m_filename.c_str(),std::ios::app);

      file << optimizer->GetCurrentIteration() << " ";
      file << optimizer->GetValue() << " ";
      file << optimizer->GetCurrentPosition()[0] << " ";
      file << optimizer->GetCurrentPosition()[1] << " ";
      file << std::endl;

      file.close();
    }

};

/** Quick function for rounding */
int roundValue(double value)
{
  return ((value < 0.0) ?
      -std::ceil( std::abs(value) -0.5 )
     : std::ceil( std::abs(value) -0.5 )
         );
}


/** Functor to compute the interferogram */
template< class TInput1, class TInput2, class TOutput>
    class SimpleInterferogramCalculator
{
  public:
  // The constructor and destructor.
    SimpleInterferogramCalculator() {};
    ~SimpleInterferogramCalculator() {};
  // Change detection operation
    inline TOutput operator()( const TInput1 & itA,
                            const TInput2 & itB)
    {

      TOutput result = itk::NumericTraits<TOutput>::Zero;
      double normA = 0.0;
      double normB = 0.0;
      for (unsigned long pos = 0; pos< itA.Size(); ++pos)
      {
        //TODO: add multiplication by G(m,i) function
        result += static_cast<TOutput>(itA.GetPixel(pos)* std::conj(itB.GetPixel(pos)));
        normA += itA.GetPixel(pos).real()*itA.GetPixel(pos).real()
            + itA.GetPixel(pos).imag()*itA.GetPixel(pos).imag();
        normB += itB.GetPixel(pos).real()*itB.GetPixel(pos).real()
            + itB.GetPixel(pos).imag()*itB.GetPixel(pos).imag();

      }
      if ((normA != 0) && (normB != 0))
        return static_cast<TOutput>( result/ (vcl_sqrt(normA)*vcl_sqrt(normB)) );
      else
        return itk::NumericTraits<TOutput>::Zero;
    }
};

int main(int argc, char* argv[])
{

  if (argc != 5)
    {
    std::cerr << "Usage: " << argv[0] << " masterImageFile slaveImageFile";
    std::cerr << " interferogram colorInterferogram" << std::endl;
    return EXIT_FAILURE;
    }

  typedef std::complex<double> PixelType;
  typedef otb::Image<PixelType,2> ImageType;
  typedef double ScalarPixelType;
  typedef otb::Image<ScalarPixelType,2> ScalarImageType;


  /* Reading master and slave images */
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

  /* Extract a region of interest */

//  unsigned int startX = 0; //ERS values
//  unsigned int startY= 2000;
  unsigned int startX = 2000; // PALSAR values
  unsigned int startY= 4000;
  unsigned int sizeX = 1000;
  unsigned int sizeY = 1000;

  typedef otb::ExtractROI<PixelType, PixelType> ExtractROIType;

  ExtractROIType::Pointer masterRoi = ExtractROIType::New();
  masterRoi->SetStartX(startX);
  masterRoi->SetStartY(startY);
  masterRoi->SetSizeX(sizeX);
  masterRoi->SetSizeY(sizeY);
  masterRoi->SetInput(master->GetOutput());
  masterRoi->UpdateOutputInformation();
  std::cout << master->GetOutput()->GetLargestPossibleRegion() << std::endl;

  ExtractROIType::Pointer slaveRoi = ExtractROIType::New();
  slaveRoi->SetStartX(startX);
//  slaveRoi->SetStartY(startY);//ERS
  slaveRoi->SetStartY(startY+376);//PALSAR
  slaveRoi->SetSizeX(sizeX);
  slaveRoi->SetSizeY(sizeY);
  slaveRoi->SetInput(slave->GetOutput());


  /* Compute Amplitude or intensity */
//  typedef otb::ComplexToIntensityImageFilter<ImageType,ScalarImageType> ComplexFilter;
  typedef itk::ComplexToModulusImageFilter<ImageType,ScalarImageType> ComplexFilter;

  ComplexFilter::Pointer masterAorI = ComplexFilter::New();
  ComplexFilter::Pointer slaveAorI = ComplexFilter::New();

  masterAorI->SetInput(masterRoi->GetOutput());
  slaveAorI->SetInput(slaveRoi->GetOutput());



  /* Smooth images to enable registration (more may be needed) */
//   typedef itk::MeanImageFilter<ScalarImageType, ScalarImageType> SmoothingFilterType;
//   typedef itk::DiscreteGaussianImageFilter<ScalarImageType, ScalarImageType> SmoothingFilterType;
  typedef otb::LeeImageFilter<ScalarImageType, ScalarImageType> SmoothingFilterType;
  SmoothingFilterType::Pointer masterFilter = SmoothingFilterType::New();
  SmoothingFilterType::Pointer slaveFilter = SmoothingFilterType::New();

  ScalarImageType::SizeType indexRadius;
  indexRadius[0] = 4; // radius along x
  indexRadius[1] = 20; // radius along y

  masterFilter->SetRadius( indexRadius );
  slaveFilter->SetRadius( indexRadius );

//   masterFilter->SetVariance(10.0);
//   slaveFilter->SetVariance(10.0);

  masterFilter->SetInput(masterAorI->GetOutput());
  slaveFilter->SetInput(slaveAorI->GetOutput());

  //Optionally output images
  typedef otb::StreamingImageFileWriter<ScalarImageType> ScalarWriterType;

  ScalarWriterType::Pointer masterWriter = ScalarWriterType:: New();
  masterWriter->SetFileName("master-filtered.tif");
  masterWriter->SetInput(masterFilter->GetOutput());
  masterWriter->Update();

  ScalarWriterType::Pointer slaveWriter = ScalarWriterType:: New();
  slaveWriter->SetFileName("slave-filtered.tif");
  slaveWriter->SetInput(slaveFilter->GetOutput());
  slaveWriter->Update();



  /* Set up registration framework */
  typedef itk::TranslationTransform< double, 2 > TransformType;
//   typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef itk::GradientDescentOptimizer       OptimizerType;

  typedef itk::MeanSquaresImageToImageMetric< ScalarImageType, ScalarImageType >    MetricType;
  typedef itk::LinearInterpolateImageFunction< ScalarImageType, double >    InterpolatorType;
  typedef itk::ImageRegistrationMethod< ScalarImageType, ScalarImageType >    RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );

  registration->SetFixedImage(    masterFilter->GetOutput()  );
  registration->SetMovingImage(   slaveFilter->GetOutput()   );


  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y

  registration->SetInitialTransformParameters( initialParameters );

//   optimizer->SetMaximumStepLength( 50 );
//   optimizer->SetMinimumStepLength( 0.01 );
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

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  CommandIterationUpdatePlot::Pointer observerPlot = CommandIterationUpdatePlot::New();
  observerPlot->SetFileName("output.dat");
  optimizer->AddObserver( itk::IterationEvent(), observerPlot );

  registration->Update();

  //second part of the registration
  registration->SetInitialTransformParameters( registration->GetLastTransformParameters());
  optimizer->SetLearningRate( 50.0 );
  optimizer->SetNumberOfIterations( 150 );
  registration->Update();

  //third part of the registration
  registration->SetInitialTransformParameters( registration->GetLastTransformParameters());
  optimizer->SetLearningRate( 30.0 );
  optimizer->SetNumberOfIterations( 50 );
  registration->Update();


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
  resample->SetInput( slaveRoi->GetOutput() );

  resample->SetSize(    masterRoi->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  masterRoi->GetOutput()->GetOrigin() );
  resample->SetOutputSpacing( masterRoi->GetOutput()->GetSpacing() );
  resample->SetDefaultPixelValue( 0 );



  typedef otb::ImageFileWriter< ImageType >  WriterFixedType;
  WriterFixedType::Pointer      writer =  WriterFixedType::New();
  writer->SetFileName( "slave-registered.tif" );
  writer->SetInput( resample->GetOutput()   );
  writer->Update();

  /* Computation of the interferogram */
  typedef SimpleInterferogramCalculator<
      itk::ConstNeighborhoodIterator<ImageType>,
      itk::ConstNeighborhoodIterator<ImageType>,
      ImageType::PixelType>  InterferogramCalculatorType;

  typedef otb::BinaryFunctorNeighborhoodImageFilter<
      ImageType,ImageType,ImageType,InterferogramCalculatorType> InterferogramFilterType;
  InterferogramFilterType::Pointer interferogram = InterferogramFilterType::New();
  interferogram->SetInput1( masterRoi->GetOutput() );
  interferogram->SetInput2( resample->GetOutput() );

  ImageType::SizeType radius;
  radius[0] = 2;
  radius[1] = 3;
  interferogram->SetRadius(radius);



  WriterType::Pointer writerInterferogram = WriterType::New();
  writerInterferogram->SetFileName(argv[3]);
  writerInterferogram->SetInput(interferogram->GetOutput());
  writerInterferogram->Update();


  /* Display the interferogram with nice colors */

  
  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef otb::Image<RGBPixelType, 2> RGBImageType;

  typedef itk::ComplexToModulusImageFilter<ImageType,ScalarImageType> ModulusFilterType;
  ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();
  modulusFilter->SetInput(masterRoi->GetOutput());

  ModulusFilterType::Pointer coherenceFilter = ModulusFilterType::New();
  coherenceFilter->SetInput(interferogram->GetOutput());

  typedef itk::ComplexToPhaseImageFilter<ImageType,ScalarImageType> PhaseFilterType;
  PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();
  phaseFilter->SetInput(interferogram->GetOutput());

  typedef otb::Functor::AmplitudePhaseToRGBFunctor
      <ScalarPixelType,ScalarPixelType,ScalarPixelType,RGBPixelType> ColorMapFunctorType;
  typedef itk::TernaryFunctorImageFilter
      <ScalarImageType, ScalarImageType, ScalarImageType, RGBImageType, ColorMapFunctorType> ColorMapFilterType;
  ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
  colormapper->GetFunctor().SetMaximum(2000);
  colormapper->GetFunctor().SetMinimum(0);


  colormapper->SetInput1(modulusFilter->GetOutput());
  colormapper->SetInput2(coherenceFilter->GetOutput());
  colormapper->SetInput3(phaseFilter->GetOutput());
  //       colormapper->SetNumberOfThreads(1);

  typedef otb::StreamingImageFileWriter<RGBImageType> WriterRGBType;
  WriterRGBType::Pointer writerRGB = WriterRGBType::New();
  writerRGB->SetFileName(argv[4]);
  writerRGB->SetInput(colormapper->GetOutput());

  writerRGB->Update();


  return EXIT_SUCCESS;
}
