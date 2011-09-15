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
#include "otbExtractROI.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "otbGenericRSTransform.h"
#include "otbLeastSquareAffineTransformEstimator.h"

#include "otbWindowedSincInterpolateImageBlackmanFunction.h"
#include "otbComplexInterpolateImageFunction.h"

#include "itkPointSet.h"
#include "otbGridIntersectionPointSetSource.h"
#include "itkFFTComplexToComplexImageFilter.h"

#include "otbStreamingResampleImageFilter.h"
#include "otbStandardWriterWatcher.h"

#include "itkFFTShiftImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkComplexToModulusImageFilter.h"

typedef std::complex< double >                                         PixelType;
typedef otb::Image< PixelType,2 >                                      ImageType;
typedef ImageType::SizeType                                            SizeType;
typedef ImageType::IndexType                                           IndexType;
typedef ImageType::RegionType                                          RegionType;

typedef itk::PointSet< PixelType, ImageType::ImageDimension >          PointSetType;
typedef otb::GridIntersectionPointSetSource< PointSetType >            PointSetSourceType;
typedef PointSetType::PointsContainer                                  PointsContainerType;
typedef PointSetType::PointType                                        PointType;

typedef otb::ImageFileReader< ImageType >                              ReaderType;
typedef otb::StreamingImageFileWriter< ImageType >                     WriterType;

typedef otb::GenericRSTransform<  >                                    TransformType; // Default is Double, 2 dimensions
typedef otb::ExtractROI< PixelType, PixelType >                        ExtractFilterType;
typedef otb::StreamingResampleImageFilter< ImageType, ImageType >      ResampleFilterType;
typedef itk::Point< PixelType::value_type,ImageType::ImageDimension >  LSQPointType;
typedef otb::LeastSquareAffineTransformEstimator< LSQPointType >       EstimateFilterType;
typedef otb::Function::BlackmanWindowFunction< PixelType::value_type > FunctionType;
typedef itk::ConstantBoundaryCondition< ImageType >                    BoundaryConditionType;
typedef PixelType::value_type                                          CoordRepType;
typedef otb::ComplexInterpolateImageFunction< ImageType,FunctionType, 
                                 BoundaryConditionType, CoordRepType > InterpolatorType;

typedef itk::FFTComplexToComplexImageFilter< PixelType::value_type, 
                                           ImageType::ImageDimension > FFTType;
typedef FFTType::OutputImageType                                       FFTOutputImageType;
typedef FFTType::TransformDirectionType                                FFTDirectionType;

typedef itk::FFTShiftImageFilter<FFTOutputImageType,FFTOutputImageType> ShiftFilterType;
typedef itk::ConstantPadImageFilter<ImageType,ImageType>          PadFilterType;
typedef otb::Image<double,2>                                           RealImageType;
typedef itk::ComplexToModulusImageFilter<FFTOutputImageType,RealImageType>      ModulusFilterType;
typedef itk::MinimumMaximumImageCalculator<RealImageType>              MinMaxCalculatorType;

//==================================== FOR VALIDATION PURPOSES ===========================================
typedef otb::ImageFileWriter<FFTOutputImageType>                       FFTWriterType;
//========================================================================================================
typedef itk::ImageRegionIteratorWithIndex< FFTOutputImageType >        ImageRegionIteratorType;

class ComplexConjugateProduct
{
public:
  inline PixelType operator()(const PixelType & a, const PixelType & b) const
  {
    return a * vcl_conj(b);
  }
};

typedef itk::BinaryFunctorImageFilter<FFTOutputImageType,FFTOutputImageType,FFTOutputImageType,ComplexConjugateProduct> ConjugateProductFilterType;

int main(int argc, char* argv[])
{

  if (argc != 12)
    {
    std::cerr << "Usage: " << argv[0] << " masterImage slaveImage tiePointsPerDim patchSizePerDim dem startx starty sizex sizey outfname1 outfname2";
    std::cerr << " OTHER PARAMS and OUTPUTS?" << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::SizeType size;
  ImageType::IndexType index;

  const char * masterImage = argv[1];
  const char * slaveImage = argv[2];
  unsigned int tiePointsPerDim = atoi(argv[3]);
  unsigned int patchSizePerDim = atoi(argv[4]);
  const char * demdir  = argv[5];
  index[0] = atoi(argv[6]);
  index[1] = atoi(argv[7]);
  size[0] = atoi(argv[8]);
  size[1] = atoi(argv[9]);
  const char * outfname1 = argv[10];
  const char * outfname2 = argv[11];

  
  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(masterImage);
  slave->SetFileName(slaveImage);

  master->GenerateOutputInformation();
  slave->GenerateOutputInformation();

  ExtractFilterType::Pointer mstExtract = ExtractFilterType::New();
  ExtractFilterType::Pointer slvExtract = ExtractFilterType::New();

  mstExtract->SetInput(master->GetOutput());
  slvExtract->SetInput(slave->GetOutput());
  
  SizeType mstSize = master->GetOutput()->GetLargestPossibleRegion().GetSize();
  
  TransformType::Pointer transform = TransformType::New();

  EstimateFilterType::Pointer estimate = EstimateFilterType::New();

  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  
  PointSetSourceType::Pointer pointSet = PointSetSourceType::New();
  PointSetType::PointType minPoint, maxPoint;
  minPoint[0] = (mstSize[0] - ((mstSize[0] - 1) % tiePointsPerDim) - 1) / tiePointsPerDim;
  minPoint[1] = (mstSize[1] - ((mstSize[1] - 1) % tiePointsPerDim) - 1) / tiePointsPerDim;
  pointSet->SetMinPoint(minPoint);
  maxPoint[0] = mstSize[0] - ((mstSize[0] - 1) % tiePointsPerDim) - 1;
  maxPoint[1] = mstSize[1] - ((mstSize[1] - 1) % tiePointsPerDim) - 1;
  pointSet->SetMaxPoint(maxPoint);

  pointSet->SetNumberOfPoints(tiePointsPerDim);
  pointSet->Update();

  // Get the the point container
  PointSetSourceType::PointsContainerPointer
  points = pointSet->GetOutput()->GetPoints();

  PointType offset;
  offset.Fill(0.0);
  ///////////////////////////////////////////////////
  // Perform genericRSTransform here if needed
  transform->SetInputKeywordList(master->GetOutput()->GetImageKeywordlist());
  transform->SetOutputKeywordList(slave->GetOutput()->GetImageKeywordlist());
  transform->SetOutputDictionary(slave->GetOutput()->GetMetaDataDictionary());
  transform->SetOutputProjectionRef(slave->GetOutput()->GetProjectionRef());
  transform->SetDEMDirectory(demdir);

  transform->InstanciateTransform();
  ///////////////////////////////////////////////////

  IndexType currentIndex;
  IndexType slvIndex;
  RegionType curMstRegion;
  RegionType curSlvRegion;
  SizeType currentSize;
  currentSize[0] = patchSizePerDim;
  currentSize[1] = patchSizePerDim;

  PointsContainerType::ConstIterator it = points->Begin();
  while (it != points->End())
  {
        PointType mstPoint = it.Value();
	PointType slvPoint = transform->TransformPoint(mstPoint);
	slvPoint[0] = floor(slvPoint[0]);
	slvPoint[1] = floor(slvPoint[1]);

	currentIndex[0] = mstPoint[0] - (patchSizePerDim / 2);
	currentIndex[1] = mstPoint[1] - (patchSizePerDim / 2);
	curMstRegion.SetIndex(currentIndex);
	curMstRegion.SetSize(currentSize);

	//slvExtract->SetExtractionRegion(currentRegion); // Should take into account GenericRSTransform-calculated initial offset (if genericRS is used)
	currentIndex[0] = slvPoint[0] - (patchSizePerDim / 2);
	currentIndex[1] = slvPoint[1] - (patchSizePerDim / 2);
	curSlvRegion.SetIndex(currentIndex);
	curSlvRegion.SetSize(currentSize);
	if(!(slave->GetOutput()->GetLargestPossibleRegion().IsInside(curSlvRegion)) || !(master->GetOutput()->GetLargestPossibleRegion().IsInside(curMstRegion)))
	{
		++it;
		continue;
	}
	mstExtract->SetExtractionRegion(curMstRegion);
	slvExtract->SetExtractionRegion(curSlvRegion);
	offset[0] = mstPoint[0] - slvPoint[0];
	offset[1] = mstPoint[1] - slvPoint[1];
	std::cout << "Initial offset: " << offset[0] << ", " << offset[1] << std::endl;

        ImageType::SizeType paddsize;
        paddsize.Fill(patchSizePerDim/2);        

        PadFilterType::Pointer pad1 = PadFilterType::New();
        pad1->SetInput(mstExtract->GetOutput());
        pad1->SetPadBound(paddsize);
        
        PadFilterType::Pointer pad2 = PadFilterType::New();
        pad2->SetInput(slvExtract->GetOutput());
        pad2->SetPadBound(paddsize);

        FFTType::Pointer fft = FFTType::New();
        fft->SetInput(pad1->GetOutput());
        fft->Update();

        FFTOutputImageType::Pointer fft1 = fft->GetOutput();

        fft = FFTType::New();
        fft->SetInput(pad2->GetOutput());
        fft->Update();

        FFTOutputImageType::Pointer fft2 = fft->GetOutput();

        ConjugateProductFilterType::Pointer conjProd = ConjugateProductFilterType::New();
        conjProd->SetInput1(fft1);
        conjProd->SetInput2(fft2);

        fft = FFTType::New();
        fft->SetInput(conjProd->GetOutput());
        fft->SetTransformDirection(FFTType::INVERSE);
        fft->Update();
  
        FFTOutputImageType::Pointer ifft = fft->GetOutput();

        ShiftFilterType::Pointer shifter = ShiftFilterType::New();
        shifter->SetInput(ifft);

        ModulusFilterType::Pointer modulus = ModulusFilterType::New();
        modulus->SetInput(shifter->GetOutput());
        modulus->Update();
  

        MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
        minMax->SetImage(modulus->GetOutput());
        minMax->ComputeMaximum();
        
        slvIndex = minMax->GetIndexOfMaximum();

	slvPoint[0] = slvPoint[0] + slvIndex[0] - (patchSizePerDim / 2);
	slvPoint[1] = slvPoint[1] + slvIndex[1] - (patchSizePerDim / 2);

	std::cout << "Master: " << mstPoint[0] << ", " << mstPoint[1];
	std::cout << " - Slave: " << slvPoint[0] << ", " << slvPoint[1] << std::endl;

	offset[0] = mstPoint[0] - slvPoint[0];
	offset[1] = mstPoint[1] - slvPoint[1];
        std::cout<<"slvIndex: "<<slvIndex<<std::endl;
	std::cout << "Final offset: " << offset[0] << ", " << offset[1] << std::endl;
    
	estimate->AddTiePoints(mstPoint, slvPoint);

        ++it;

	// if(mstPoint[0] == 4364 && mstPoint[1] == 2472)
	// {
        // WriterType::Pointer mstWriter = WriterType::New();
        // mstWriter->SetInput(mstExtract->GetOutput());
        // mstWriter->SetFileName("mst.tif");
        // mstWriter->Update();

        // WriterType::Pointer slvWriter = WriterType::New();
        // slvWriter->SetInput(slvExtract->GetOutput());
        // slvWriter->SetFileName("slv.tif");
        // slvWriter->Update();

        // FFTWriterType::Pointer mstFFTWriter = FFTWriterType::New();
        // mstFFTWriter->SetInput(fft1);
        // mstFFTWriter->SetFileName("masterFFT.tif");
        // mstFFTWriter->Update();
        
        // FFTWriterType::Pointer slvFFTWriter = FFTWriterType::New();
        // slvFFTWriter->SetInput(fft2);
        // slvFFTWriter->SetFileName("slaveFFT.tif");
        // slvFFTWriter->Update();
        
        // FFTWriterType::Pointer crossFFTWriter = FFTWriterType::New();
        // crossFFTWriter->SetInput(shifter->GetOutput());
        // crossFFTWriter->SetFileName("invcrossFFT.tif");
        // crossFFTWriter->Update();
	// }
  }

  estimate->Compute();

  EstimateFilterType::CovariantVectorType rmsError = estimate->GetRMSError();
  EstimateFilterType::CovariantVectorType relResidual = estimate->GetRelativeResidual();

  std::cout << "RMS error is:" << rmsError[0] << " in range and " << rmsError[1] << " in azimuth." << std::endl;
  std::cout << "Relative residual is:" << relResidual[0] << " in range and " << relResidual[1] << " in azimuth." << std::endl;
  
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(slave->GetOutput());
  interpolator->SetRadius(3);
  interpolator->SetNormalizeZeroFrequency(0.01);

  resample->SetTransform(estimate->GetAffineTransform());
  resample->SetInterpolator(interpolator);
  resample->SetInput(slave->GetOutput());
  resample->SetOutputSize(mstSize);
  resample->SetOutputOrigin(master->GetOutput()->GetOrigin());
  resample->SetOutputSpacing(master->GetOutput()->GetSpacing());


  ExtractFilterType::Pointer extract = ExtractFilterType::New();
  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(index);
  extract->SetExtractionRegion(region);
  extract->SetInput(master->GetOutput());

  typedef otb::StreamingImageFileWriter< ImageType > WriterFixedType;
  WriterFixedType::Pointer writer = WriterFixedType::New();
  writer->SetFileName(outfname1);
  writer->SetInput(extract->GetOutput());

  otb::StandardWriterWatcher watcher1(writer,resample,"Extracting from master image.");

  writer->Update();

  extract = ExtractFilterType::New();
  extract->SetExtractionRegion(region);
  extract->SetInput(resample->GetOutput());
  
  writer = WriterFixedType::New();
  writer->SetFileName(outfname2);
  writer->SetInput(extract->GetOutput());

  otb::StandardWriterWatcher watcher2(writer,resample,"Extracting from registered slave image.");
  writer->Update();

  return EXIT_SUCCESS;
}

