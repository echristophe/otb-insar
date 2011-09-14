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

//==================================== FOR VALIDATION PURPOSES ===========================================
typedef otb::ImageFileWriter<FFTOutputImageType>                       FFTWriterType;
//========================================================================================================
typedef itk::ImageRegionIteratorWithIndex< FFTOutputImageType >        ImageRegionIteratorType;

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

  FFTOutputImageType::Pointer crossImage = FFTOutputImageType::New();

  currentIndex[0] = 0;
  currentIndex[1] = 0;
  curMstRegion.SetIndex(currentIndex);
  curMstRegion.SetSize(currentSize);
  crossImage->SetSpacing(master->GetOutput()->GetSpacing());
  crossImage->SetOrigin(master->GetOutput()->GetOrigin());
  crossImage->SetRegions(curMstRegion);
  crossImage->Allocate();

  FFTType::Pointer crossFFT = FFTType::New();
  crossFFT->SetTransformDirection(FFTType::INVERSE);

  PointsContainerType::ConstIterator it = points->Begin();
  while (it != points->End())
  {
	crossImage->FillBuffer(0.0);

        PointType mstPoint = it.Value();
	PointType slvPoint = transform->TransformPoint(mstPoint);
	slvPoint[0] = floor(slvPoint[0]);
	slvPoint[1] = floor(slvPoint[1]);
	//if(mstPoint[1] > 5000)
		//int j = 1;
		//	break;
	currentIndex[0] = mstPoint[0] - (patchSizePerDim / 2);
	currentIndex[1] = mstPoint[1] - (patchSizePerDim / 2);
	curMstRegion.SetIndex(currentIndex);
	curMstRegion.SetSize(currentSize);
	//currentRegion.Crop(master->GetOutput()->GetLargestPossibleRegion());

	//slvExtract->SetExtractionRegion(currentRegion); // Should take into account GenericRSTransform-calculated initial offset (if genericRS is used)
	currentIndex[0] = slvPoint[0] - (patchSizePerDim / 2);
	currentIndex[1] = slvPoint[1] - (patchSizePerDim / 2);
	curSlvRegion.SetIndex(currentIndex);
	curSlvRegion.SetSize(currentSize);
	//currentRegion.Crop(slave->GetOutput()->GetLargestPossibleRegion());
	//if(!currentRegion.IsInside(slave->GetOutput()->GetLargestPossibleRegion()))
	if(!(slave->GetOutput()->GetLargestPossibleRegion().IsInside(curSlvRegion)) || !(master->GetOutput()->GetLargestPossibleRegion().IsInside(curMstRegion)))
	{
		/*
		currentIndex[0] = mstPoint[0];
		currentIndex[1] = mstPoint[1];
		currentRegion.SetIndex(currentIndex);
		currentRegion.SetSize(currentSize);
		currentRegion.Crop(slave->GetOutput()->GetLargestPossibleRegion());
		*/
		++it;
		continue;
	}
	mstExtract->SetExtractionRegion(curMstRegion);
	slvExtract->SetExtractionRegion(curSlvRegion);
	offset[0] = mstPoint[0] - slvPoint[0];
	offset[1] = mstPoint[1] - slvPoint[1];
	std::cout << "Initial offset: " << offset[0] << ", " << offset[1] << std::endl;

	//mstFFT->Update();
	//slvFFT->Update();
	mstFFT->UpdateLargestPossibleRegion();
	slvFFT->UpdateLargestPossibleRegion();

	ImageRegionIteratorType mstIt(mstFFT->GetOutput(), mstFFT->GetOutput()->GetRequestedRegion());
	ImageRegionIteratorType slvIt(slvFFT->GetOutput(), slvFFT->GetOutput()->GetRequestedRegion());

	ImageRegionIteratorType crossIt(crossImage, crossImage->GetRequestedRegion());

	for(mstIt.GoToBegin(),	slvIt.GoToBegin(), crossIt.GoToBegin(); !mstIt.IsAtEnd(), !slvIt.IsAtEnd(); ++mstIt, ++slvIt, ++crossIt)
	{
		crossIt.Value() = mstIt.Value()*conj(slvIt.Value());
	}

	crossFFT->SetInput(crossImage);
	crossFFT->Update();

	ImageRegionIteratorType invIt(crossFFT->GetOutput(), crossFFT->GetOutput()->GetRequestedRegion());

	double maxValue = 0.0;
	//PointType slvPoint = mstPoint;

	slvIndex.Fill(0.0);
	for(invIt.GoToBegin(); !invIt.IsAtEnd(); ++invIt)
	{
		double value = abs(invIt.Value());

		if(value > maxValue)
		{
			maxValue = value;

			slvIndex = invIt.GetIndex();

			//slvPoint[0] = slvIndex[0] + mstPoint[0] + offset[0];
			//slvPoint[1] = slvIndex[1] + mstPoint[1] + offset[1];			
		}
	}
	slvPoint[0] = slvPoint[0] + slvIndex[0] - (patchSizePerDim / 2);
	slvPoint[1] = slvPoint[1] + slvIndex[1] - (patchSizePerDim / 2);

	std::cout << "Master: " << mstPoint[0] << ", " << mstPoint[1];
	std::cout << " - Slave: " << slvPoint[0] << ", " << slvPoint[1] << std::endl;

	offset[0] = mstPoint[0] - slvPoint[0];
	offset[1] = mstPoint[1] - slvPoint[1];
	std::cout << "Final offset: " << offset[0] << ", " << offset[1] << std::endl;
    
	estimate->AddTiePoints(mstPoint, slvPoint);

    ++it;

	// if(mstPoint[0] == 5100 && mstPoint[1] == 3474)
	// {
	// 	FFTWriterType::Pointer mstFFTWriter = FFTWriterType::New();
	// 	mstFFTWriter->SetInput(mstFFT->GetOutput());
	// 	mstFFTWriter->SetFileName("masterFFTx5100_y3474_256.hdr");
	// 	mstFFTWriter->Update();
	// 	FFTWriterType::Pointer slvFFTWriter = FFTWriterType::New();
	// 	slvFFTWriter->SetInput(slvFFT->GetOutput());
	// 	slvFFTWriter->SetFileName("slaveFFTx5100_y3474_256.hdr");
	// 	slvFFTWriter->Update();
	// 	FFTWriterType::Pointer crossFFTWriter = FFTWriterType::New();
	// 	crossFFTWriter->SetInput(crossFFT->GetOutput());
	// 	crossFFTWriter->SetFileName("crossFFTx5100_y3474_256.hdr");
	// 	crossFFTWriter->Update();
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

