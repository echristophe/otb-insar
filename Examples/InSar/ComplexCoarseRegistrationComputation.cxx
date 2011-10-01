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
#include "itkDivideImageFilter.h"

#include "itkConstNeighborhoodIterator.h"

#include "otbSubsampledImageRegionConstIterator.h"
#include "itkAffineTransform.h"
#include "otbStreamingWarpImageFilter.h"

#include "itkLinearInterpolateImageFunction.h"

typedef std::complex< double >                                         PixelType;
typedef otb::Image< PixelType,2 >                                      ImageType;
typedef ImageType::SizeType                                            SizeType;
typedef ImageType::IndexType                                           IndexType;
typedef ImageType::RegionType                                          RegionType;
typedef ImageType::SpacingType										   SpacingType;

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
typedef InterpolatorType::ContinuousIndexType            ContinuousIndexType;

typedef itk::FFTComplexToComplexImageFilter< PixelType::value_type, 
                                           ImageType::ImageDimension > FFTType;
typedef FFTType::OutputImageType                                       FFTOutputImageType;
typedef FFTType::TransformDirectionType                                FFTDirectionType;

typedef itk::FFTShiftImageFilter<FFTOutputImageType,
												   FFTOutputImageType> ShiftFilterType;
typedef itk::ConstantPadImageFilter<ImageType,ImageType>			   PadFilterType;
typedef otb::Image<double,2>										   RealImageType;
typedef itk::ComplexToModulusImageFilter<FFTOutputImageType,
													RealImageType>	   ModulusFilterType;
typedef itk::DivideImageFilter<FFTOutputImageType,RealImageType,
												FFTOutputImageType>    DivideFilterType;
typedef itk::MinimumMaximumImageCalculator<RealImageType>			   MinMaxCalculatorType;

//==================================== FOR VALIDATION PURPOSES ===========================================
typedef otb::ImageFileWriter<FFTOutputImageType>                       FFTWriterType;
//========================================================================================================
typedef itk::ImageRegionIteratorWithIndex< FFTOutputImageType >        ImageRegionIteratorType;
typedef itk::ConstNeighborhoodIterator< ImageType >					   ImageConstNeighborhoodIteratorType;

typedef otb::SubsampledImageRegionConstIterator< ImageType >		   SubsampledImageRegionConstIteratorType;
typedef itk::Vector<double, 2>										   DeformationOffsetPixelType;
typedef otb::Image<DeformationOffsetPixelType, 2>                      DeformationFieldType;
typedef DeformationFieldType::RegionType							   DeformationRegionType;
typedef itk::ImageRegionIterator< DeformationFieldType >			   DeformationImageRegionIteratorType;
typedef otb::StreamingImageFileWriter< DeformationFieldType >		   DeformationImageWriterType;

typedef EstimateFilterType::AffineTransformType						   AffineTransformType;
typedef otb::StreamingWarpImageFilter< ImageType, ImageType, 
	DeformationFieldType>											   WarpImageFilterType;

typedef otb::StreamingResampleImageFilter< DeformationFieldType, 
											DeformationFieldType >     DeformationResampleFilterType;
typedef itk::LinearInterpolateImageFunction< DeformationFieldType >    LinearInterpolateFilterType;

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

  if (argc != 18)
    {
    std::cerr << "Usage: " << argv[0] << " masterImage slaveImage tiePointsPerDim patchSizePerDim dem startx starty sizex sizey outfname1 outfname2 coherency_threshold subPixelAccuracy subSampleFactor windowSize fine_threshold outfname3";
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
  double coherency_threshold = atof(argv[12]);
  double subPixelAccuracy = atof(argv[13]);
  IndexType::IndexValueType subSampleFactor = atoi(argv[14]);
  unsigned int windowSizePerDim = atoi(argv[15]);
  double fine_threshold = atof(argv[16]);
  const char * outfname3 = argv[17];

  
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

  // Prepare Generic RS Transform
  transform->SetInputKeywordList(master->GetOutput()->GetImageKeywordlist());
  transform->SetOutputKeywordList(slave->GetOutput()->GetImageKeywordlist());
  transform->SetDEMDirectory(demdir);

  transform->InstanciateTransform();
  

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

        // Try to normalise to get the normalized coherence coeff
        ModulusFilterType::Pointer conjProdMod = ModulusFilterType::New();
        conjProdMod->SetInput(conjProd->GetOutput());

        DivideFilterType::Pointer conjProdNorm = DivideFilterType::New();
        conjProdNorm->SetInput1(conjProd->GetOutput());
        conjProdNorm->SetInput2(conjProdMod->GetOutput());
        
        fft = FFTType::New();
        fft->SetInput(conjProdNorm->GetOutput());
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

	slvPoint[0] = slvPoint[0] - slvIndex[0] + (patchSizePerDim / 2);
	slvPoint[1] = slvPoint[1] - slvIndex[1] + (patchSizePerDim / 2);

	std::cout << "Master: " << mstPoint[0] << ", " << mstPoint[1];
	std::cout << " - Slave: " << slvPoint[0] << ", " << slvPoint[1] << std::endl;
	std::cout << "Final offset: " << slvIndex[0] - (patchSizePerDim / 2.) << ", " << slvIndex[1] - (patchSizePerDim / 2.) << std::endl;
        std::cout<<"Coherency: "<<minMax->GetMaximum()<<std::endl;
        
        if(minMax->GetMaximum()>coherency_threshold)
          {
          estimate->AddTiePoints(mstPoint, slvPoint);
          std::cout<<"Tie points added"<<std::endl;
          }
        std::cout<<"==================================="<<std::endl;
        ++it;
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

  std::cout<<"Transform matrix: "<<std::endl;
  std::cout<<estimate->GetAffineTransform()->GetMatrix()<<std::endl;

  std::cout<<"Transform offset: "<<std::endl;
  std::cout<<estimate->GetAffineTransform()->GetTranslation()<<std::endl;

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

  otb::StandardWriterWatcher watcher1(writer,extract,"Extracting from master image.");

  writer->Update();

  extract = ExtractFilterType::New();
  extract->SetExtractionRegion(region);
  extract->SetInput(resample->GetOutput());
  
  writer = WriterFixedType::New();
  writer->SetFileName(outfname2);
  writer->SetInput(extract->GetOutput());

  otb::StandardWriterWatcher watcher2(writer,resample,"Extracting from registered slave image.");
  writer->Update();

  /** =================================================================================== */
  /** Fine registration starts here ===================================================== */
  //Get grid iterator for master image
  DeformationRegionType::IndexType defIndex;
  defIndex.Fill(0);
  DeformationRegionType::SizeType defSize;
  defSize[0] = (master->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - (master->GetOutput()->GetLargestPossibleRegion().GetSize()[0] % subSampleFactor)) / subSampleFactor;
  defSize[1] = (master->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - (master->GetOutput()->GetLargestPossibleRegion().GetSize()[1] % subSampleFactor)) / subSampleFactor;

  DeformationRegionType deformationRegion;
  deformationRegion.SetIndex(defIndex);
  deformationRegion.SetSize(defSize);

  DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
  deformationField->SetRegions(deformationRegion);
  deformationField->SetSpacing(master->GetOutput()->GetSpacing() * subSampleFactor);
  deformationField->Allocate();

  DeformationImageRegionIteratorType deformationIt(deformationField, deformationField->GetRequestedRegion());

  AffineTransformType::Pointer affineTransformFilter = estimate->GetAffineTransform();

  currentSize.Fill(windowSizePerDim);

  IndexType mstIndex;
  mstIndex.Fill(0);
  //LOOP
  deformationIt.GoToBegin();
  while(!deformationIt.IsAtEnd())
  {
	//Get physical point from index

	  PointType masterPoint;
	  master->GetOutput()->TransformIndexToPhysicalPoint(mstIndex, masterPoint);

	//Transform physical point in master to physical point in slave
	  PointType slavePoint = affineTransformFilter->TransformPoint(masterPoint);

	//Get slave index from physical point

	  IndexType slvIndex;
	  slave->GetOutput()->TransformPhysicalPointToIndex(slavePoint, slvIndex);

	  DeformationOffsetPixelType deformationOffset;

	  deformationOffset[0] = slavePoint[0] - masterPoint[0];
	  deformationOffset[1] = slavePoint[1] - masterPoint[1];

	//Get region around point in master and slave

	  IndexType tmpMstIndex;
	  tmpMstIndex[0] = mstIndex[0] - (windowSizePerDim / 2);
	  tmpMstIndex[1] = mstIndex[1] - (windowSizePerDim / 2);
	  curMstRegion.SetIndex(tmpMstIndex);
	  curMstRegion.SetSize(currentSize);

	  slvIndex[0] = slvIndex[0] - (windowSizePerDim / 2);
	  slvIndex[1] = slvIndex[1] - (windowSizePerDim / 2);
	  curSlvRegion.SetIndex(slvIndex);
	  curSlvRegion.SetSize(currentSize);

	  if(!(slave->GetOutput()->GetLargestPossibleRegion().IsInside(curSlvRegion)) || !(master->GetOutput()->GetLargestPossibleRegion().IsInside(curMstRegion)))
	  {
		  deformationIt.Set(deformationOffset);

		  mstIndex[0] = mstIndex[0] + subSampleFactor;
		  if(mstIndex[0] > master->GetOutput()->GetLargestPossibleRegion().GetSize()[0])
		  {
			  mstIndex[0] = 0;
			  mstIndex[1] = mstIndex[1] + subSampleFactor;
		  }

		  ++deformationIt;
		  continue;
	  }
	  mstExtract->SetExtractionRegion(curMstRegion);
	  slvExtract->SetExtractionRegion(curSlvRegion);

	//Pad region to x times original size to get 1/x pixel precision

	  ImageType::SizeType paddsize;
      paddsize.Fill(((windowSizePerDim/subPixelAccuracy) - windowSizePerDim) / 2);        

      PadFilterType::Pointer pad1 = PadFilterType::New();
      pad1->SetInput(mstExtract->GetOutput());
      pad1->SetPadBound(paddsize);
        
      PadFilterType::Pointer pad2 = PadFilterType::New();
      pad2->SetInput(slvExtract->GetOutput());
      pad2->SetPadBound(paddsize);

	//Perform fft on both image regions

	  FFTType::Pointer fft = FFTType::New();
      fft->SetInput(pad1->GetOutput());
      fft->Update();

      FFTOutputImageType::Pointer fft1 = fft->GetOutput();

      fft = FFTType::New();
      fft->SetInput(pad2->GetOutput());
      fft->Update();

      FFTOutputImageType::Pointer fft2 = fft->GetOutput();

	//Calculate correlation (complex conjugate product)
	
	  ConjugateProductFilterType::Pointer conjProd = ConjugateProductFilterType::New();
      conjProd->SetInput1(fft1);
      conjProd->SetInput2(fft2);

	//Normalize

	  ModulusFilterType::Pointer conjProdMod = ModulusFilterType::New();
      conjProdMod->SetInput(conjProd->GetOutput());

      DivideFilterType::Pointer conjProdNorm = DivideFilterType::New();
      conjProdNorm->SetInput1(conjProd->GetOutput());
      conjProdNorm->SetInput2(conjProdMod->GetOutput());

	//Get maximum value index

	  fft = FFTType::New();
      fft->SetInput(conjProdNorm->GetOutput());
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
      
	  if(minMax->GetMaximum() > fine_threshold)
	  {
		IndexType tmpIndex = minMax->GetIndexOfMaximum();

		//Transform index to physical tie point

        deformationOffset[0] = deformationOffset[0] + (double(tmpIndex[0]) - windowSizePerDim / 2) * subPixelAccuracy;
		deformationOffset[1] = deformationOffset[1] + (double(tmpIndex[1]) - windowSizePerDim / 2) * subPixelAccuracy;
	  }

	//Add tie point to displacement map

	  deformationIt.Set(deformationOffset);

	  mstIndex[0] = mstIndex[0] + subSampleFactor;
	  if(mstIndex[0] > master->GetOutput()->GetLargestPossibleRegion().GetSize()[0])
	  {
		  mstIndex[0] = 0;
		  mstIndex[1] = mstIndex[1] + subSampleFactor;
	  }

	  std::cout << "Master: " << masterPoint[0] << ", " << masterPoint[1];
	  std::cout << " - Slave: " << slavePoint[0] << ", " << slavePoint[1] << std::endl;
	  std::cout << "Final offset: " << deformationOffset[0] << ", " << deformationOffset[1] << std::endl;
      std::cout<<"Coherency: "<<minMax->GetMaximum()<<std::endl;

	  ++deformationIt;
  }

  DeformationResampleFilterType::Pointer resampleDefField = DeformationResampleFilterType::New();

  resampleDefField->SetInput(deformationField);

  LinearInterpolateFilterType::Pointer linearInterpolator = LinearInterpolateFilterType::New();

  resampleDefField->SetInterpolator(linearInterpolator);
  resampleDefField->SetOutputSize(slave->GetOutput()->GetLargestPossibleRegion().GetSize());
  resampleDefField->SetOutputOrigin(slave->GetOutput()->GetOrigin());
  resampleDefField->SetOutputSpacing(slave->GetOutput()->GetSpacing());
/*
  DeformationImageWriterType::Pointer writeDefField = DeformationImageWriterType::New();
  writeDefField->SetInput(resampleDefField->GetOutput());
  writeDefField->SetFileName("fineDefField.hdr");
  writeDefField->SetNumberOfDivisionsTiledStreaming(1);
  writeDefField->Update();
*/
  WarpImageFilterType::Pointer warp = WarpImageFilterType::New();

  warp->SetInput(slave->GetOutput());
  //warp->SetDeformationField(deformationField);  // Not currently working if grid step is different from 1
  warp->SetDeformationField(resampleDefField->GetOutput());
  warp->SetEdgePaddingValue(0.0);

  DeformationOffsetPixelType maxDeformation;
  maxDeformation[0] = estimate->GetAffineTransform()->GetTranslation()[0] - 10.0; // The minus is offset related... need a if statement!
  maxDeformation[1] = estimate->GetAffineTransform()->GetTranslation()[1] - 10.0; // The minus is offset related... same as above!

  warp->SetMaximumDeformation(maxDeformation);
  warp->SetInterpolator(interpolator);

  extract = ExtractFilterType::New();
  extract->SetExtractionRegion(region);
  extract->SetInput(warp->GetOutput());
  
  writer = WriterFixedType::New();
  writer->SetFileName(outfname3);
  writer->SetInput(extract->GetOutput());
  //writer->SetInput(warp->GetOutput());

  otb::StandardWriterWatcher watcher3(writer,warp,"Extracting from fine registered slave image.");

  writer->Update();


	/** Fine registration ends here ======================================================= */
	/** =================================================================================== */

  return EXIT_SUCCESS;
}

