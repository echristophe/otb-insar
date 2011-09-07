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

#include "otbComplexInterpolateImageFunction.h"
#include "otbWindowedSincInterpolateImageBlackmanFunction.h"

#include "itkPointSet.h"
#include "otbGridIntersectionPointSetSource.h"
#include "itkFFTComplexToComplexImageFilter.h"
#include "otbStreamingResampleImageFilter.h"

#include "otbStandardWriterWatcher.h"

typedef std::complex<double>						PixelType;
typedef otb::Image<PixelType,2>						ImageType;
typedef ImageType::SizeType							SizeType;
typedef ImageType::IndexType						IndexType;
typedef ImageType::RegionType						RegionType;

typedef otb::ImageFileReader<ImageType>				ReaderType;
typedef otb::StreamingImageFileWriter<ImageType>	WriterType;

typedef otb::GenericRSTransform< >					TransformType; // Default is Double, 2 dimensions

typedef itk::PointSet<PixelType, 2>							PointSetType;
typedef otb::GridIntersectionPointSetSource<PointSetType>	PointSetSourceType;
typedef PointSetType::PointsContainer						PointsContainerType;
typedef PointSetType::PointType								PointType;

typedef otb::ExtractROI< PixelType, PixelType >				ExtractFilterType;

typedef itk::FFTComplexToComplexImageFilter<double, 2>			FFTType;
typedef FFTType::OutputImageType								FFTOutputImageType;

typedef itk::ImageRegionIteratorWithIndex<FFTOutputImageType>	ImageRegionIteratorType;

typedef otb::StreamingResampleImageFilter< ImageType, ImageType >		ResampleFilterType;

typedef itk::Point<double,2>						myPointType;
typedef otb::LeastSquareAffineTransformEstimator<myPointType>		EstimateFilterType;

typedef otb::Function::BlackmanWindowFunction<double> FunctionType;
typedef itk::ConstantBoundaryCondition<ImageType> BoundaryConditionType;
typedef double CoordRepType;
typedef otb::ComplexInterpolateImageFunction<ImageType,FunctionType, BoundaryConditionType, CoordRepType> InterpolatorType;



int main(int argc, char* argv[])
{

  if (argc != 6)
    {
    std::cerr << "Usage: " << argv[0] << " masterImage slaveImage tiePointsPerDim patchSizePerDim registeredSlaveImageFile";
    std::cerr << " OTHER PARAMS and OUTPUTS?" << std::endl;
    return EXIT_FAILURE;
    }

  const char * masterImage = argv[1];
  const char * slaveImage = argv[2];
  unsigned int tiePointsPerDim = atoi(argv[3]);
  unsigned int patchSizePerDim = atoi(argv[4]);
  const char * outRegisteredSlave = argv[5];
  
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

  FFTType::Pointer mstFFT = FFTType::New();
  FFTType::Pointer slvFFT = FFTType::New();

  mstFFT->SetInput(mstExtract->GetOutput());
  slvFFT->SetInput(slvExtract->GetOutput());
  
  SizeType mstSize = master->GetOutput()->GetLargestPossibleRegion().GetSize();
  
  TransformType::Pointer transform = TransformType::New();

  EstimateFilterType::Pointer estimate = EstimateFilterType::New();

  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  
  PointSetSourceType::Pointer pointSet = PointSetSourceType::New();
  PointSetType::PointType minPoint, maxPoint;
  minPoint[0] = (mstSize[0] - (mstSize[0] % tiePointsPerDim)) / tiePointsPerDim;
  minPoint[1] = (mstSize[1] - (mstSize[1] % tiePointsPerDim)) / tiePointsPerDim;
  pointSet->SetMinPoint(minPoint);
  maxPoint[0] = mstSize[0] - (mstSize[0] % tiePointsPerDim);
  maxPoint[1] = mstSize[1] - (mstSize[1] % tiePointsPerDim);
  pointSet->SetMaxPoint(maxPoint);

  pointSet->SetNumberOfPoints(tiePointsPerDim);
  pointSet->Update();

  // Get the the point container
  PointSetSourceType::PointsContainerPointer
  points = pointSet->GetOutput()->GetPoints();
/*  
  // Temp test to insure no "in between points" at non existing indexes
  PointsContainerType::ConstIterator it = points->Begin();
  while (it != points->End())
    {
    PointSetType::PointType p = it.Value();
    std::cout.width(5); std::cout << p[0] << ", ";
    std::cout.width(5); std::cout << p[1] << std::endl;
    ++it;
    }
*/

  PointType offset;
  offset.Fill(0.0);
  ///////////////////////////////////////////////////
  // Perform genericRSTransform here if needed
  ///////////////////////////////////////////////////

  IndexType currentIndex;
  IndexType slvIndex;
  RegionType currentRegion;
  SizeType currentSize;
  currentSize[0] = patchSizePerDim;
  currentSize[1] = patchSizePerDim;

  PointsContainerType::ConstIterator it = points->Begin();
  while (it != points->End())
  {
    PointType mstPoint = it.Value();
	if(mstPoint[1] > 5983)
		break;
	currentIndex[0] = mstPoint[0];
	currentIndex[1] = mstPoint[1];
	currentRegion.SetIndex(currentIndex);
	currentRegion.SetSize(currentSize);
	currentRegion.Crop(master->GetOutput()->GetLargestPossibleRegion());

	mstExtract->SetExtractionRegion(currentRegion);
	slvExtract->SetExtractionRegion(currentRegion); // Should take into account GenericRSTransform-calculated initial offset (if genericRS is used)
/*
	mstExtract->SetStartX(mstPoint[0]);
	mstExtract->SetStartY(mstPoint[1]);
	mstExtract->SetSizeX(patchSizePerDim);
	mstExtract->SetSizeY(patchSizePerDim);

	slvExtract->SetStartX(mstPoint[0] + offset[0]);
	slvExtract->SetStartY(mstPoint[1] + offset[1]);
	slvExtract->SetSizeX(patchSizePerDim);
	slvExtract->SetSizeY(patchSizePerDim);
*/
	mstFFT->Update();
	slvFFT->Update();

	ImageRegionIteratorType mstIt(mstFFT->GetOutput(), mstFFT->GetOutput()->GetRequestedRegion());
	ImageRegionIteratorType slvIt(slvFFT->GetOutput(), slvFFT->GetOutput()->GetRequestedRegion());

	double maxValue = 0.0;
	PointType slvPoint = mstPoint;

	for(mstIt.GoToBegin(),	slvIt.GoToBegin(); !mstIt.IsAtEnd(), !slvIt.IsAtEnd(); ++mstIt, ++slvIt)
	{
		double value = abs(mstIt.Value()*conj(slvIt.Value()));

		if(value > maxValue)
		{
			maxValue = value;

			slvIndex = slvIt.GetIndex();

			slvPoint[0] = slvIndex[0] + mstPoint[0] + offset[0];
			slvPoint[1] = slvIndex[1] + mstPoint[1] + offset[1];
		}
	}
	estimate->AddTiePoints(mstPoint, slvPoint);

    ++it;
  }

  estimate->Compute();

  EstimateFilterType::CovariantVectorType rmsError = estimate->GetRMSError();
  EstimateFilterType::CovariantVectorType relResidual = estimate->GetRelativeResidual();

  std::cout << "RMS error is:" << rmsError[0] << " in range and " << rmsError[1] << " in azimuth." << std::endl;
  std::cout << "Relative residual is:" << relResidual[0] << " in range and " << relResidual[1] << " in azimuth." << std::endl;

    // Set-up interpolator
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(slave->GetOutput());
  interpolator->SetRadius(3);
  interpolator->SetNormalizeZeroFrequency(0.01);
  resample->SetInterpolator(interpolator);
  
  resample->SetTransform(estimate->GetAffineTransform());
  resample->SetInput(slave->GetOutput());
  resample->SetOutputSize(mstSize);
  resample->SetOutputOrigin(master->GetOutput()->GetOrigin());
  resample->SetOutputSpacing(master->GetOutput()->GetSpacing());

  typedef otb::StreamingImageFileWriter<ImageType> WriterFixedType;
  WriterFixedType::Pointer writer = WriterFixedType::New();
  writer->SetFileName(outRegisteredSlave);
  writer->SetInput(resample->GetOutput());

  otb::StandardWriterWatcher watcher1(writer,resample,"Resampling slave image.");

  writer->Update();

  return EXIT_SUCCESS;
}

