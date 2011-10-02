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

#include "itkFFTComplexToComplexImageFilter.h"

#include "otbStandardWriterWatcher.h"

#include "itkFFTShiftImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"

#include "itkMeanImageFunction.h"

typedef std::complex< double >                                         PixelType;
typedef otb::Image< PixelType,2 >                                      ImageType;
typedef ImageType::SizeType                                            SizeType;
typedef ImageType::IndexType                                           IndexType;
typedef ImageType::RegionType                                          RegionType;
typedef ImageType::SpacingType										   SpacingType;

typedef otb::ImageFileReader< ImageType >                              ReaderType;
typedef otb::StreamingImageFileWriter< ImageType >                     WriterType;

typedef otb::ExtractROI< PixelType, PixelType >                        ExtractFilterType;

typedef itk::FFTComplexToComplexImageFilter< PixelType::value_type, 
                                           ImageType::ImageDimension > FFTType;
typedef FFTType::OutputImageType                                       FFTOutputImageType;

typedef itk::FFTShiftImageFilter<FFTOutputImageType,FFTOutputImageType> ShiftFilterType;
typedef itk::ConstantPadImageFilter<ImageType,ImageType>          PadFilterType;
typedef otb::Image<double,2>                                           RealImageType;
typedef itk::ComplexToModulusImageFilter<FFTOutputImageType,RealImageType>      ModulusFilterType;
typedef itk::DivideImageFilter<FFTOutputImageType,RealImageType,FFTOutputImageType> DivideFilterType;
typedef itk::MinimumMaximumImageCalculator<RealImageType>              MinMaxCalculatorType;

typedef itk::ImageRegionIteratorWithIndex< ImageType >					ImageRegionIteratorType;
typedef itk::ImageRegionIteratorWithIndex< RealImageType >				RealImageRegionIteratorType;
typedef itk::ConstNeighborhoodIterator< ImageType >					   ImageConstNeighborhoodIteratorType;

typedef itk::MeanImageFunction< RealImageType, PixelType::value_type > MeanImageFunctionType;

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

  if (argc != 6)
    {
    std::cerr << "Usage: " << argv[0] << " masterImage slaveImage patchSizeRange patchSizeAzimuth outfname1";
    return EXIT_FAILURE;
    }
  ImageType::SizeType size;
  ImageType::IndexType index;

  const char * masterImage = argv[1];
  const char * slaveImage = argv[2];
  unsigned int patchSizeRange = atoi(argv[3]);
  unsigned int patchSizeAzimuth = atoi(argv[4]);
  const char * outfname1 = argv[5];
    
  ReaderType::Pointer master = ReaderType::New();
  ReaderType::Pointer slave = ReaderType::New();

  master->SetFileName(masterImage);
  slave->SetFileName(slaveImage);

  master->GenerateOutputInformation();
  slave->GenerateOutputInformation();
  
  SizeType mstSize = master->GetOutput()->GetLargestPossibleRegion().GetSize();

  RegionType curMstRegion;
  RegionType curSlvRegion;
  SizeType regionSize;
  regionSize[0] = patchSizeRange;
  regionSize[1] = patchSizeAzimuth;

  RealImageType::Pointer coherenceImage = RealImageType::New();
  coherenceImage->SetRegions(master->GetOutput()->GetLargestPossibleRegion());
  coherenceImage->SetSpacing(master->GetOutput()->GetSpacing());
  coherenceImage->SetOrigin(master->GetOutput()->GetOrigin());
  coherenceImage->Allocate();

  ExtractFilterType::Pointer mstExtract = ExtractFilterType::New();
  mstExtract->SetInput(master->GetOutput());

  ExtractFilterType::Pointer slvExtract = ExtractFilterType::New();
  slvExtract->SetInput(slave->GetOutput());

  MeanImageFunctionType::Pointer mean = MeanImageFunctionType::New();

  mean->SetNeighborhoodRadius(1.0);

  master->Update();

  ImageRegionIteratorType mstIt(master->GetOutput(), master->GetOutput()->GetLargestPossibleRegion());
  RealImageRegionIteratorType coIt(coherenceImage, coherenceImage->GetLargestPossibleRegion());
  for(mstIt.GoToBegin(), coIt.GoToBegin(); !mstIt.IsAtEnd(); ++mstIt, ++coIt)
  {
    IndexType mstIndex = mstIt.GetIndex();

	mstIndex[0] = mstIndex[0] - (patchSizeRange / 2);
	mstIndex[1] = mstIndex[1] - (patchSizeAzimuth / 2);
	curMstRegion.SetIndex(mstIndex);
	curMstRegion.SetSize(regionSize);

	curSlvRegion.SetIndex(mstIndex);
	curSlvRegion.SetSize(regionSize);

	if(!(slave->GetOutput()->GetLargestPossibleRegion().IsInside(curSlvRegion)) || !(master->GetOutput()->GetLargestPossibleRegion().IsInside(curMstRegion)))
	{
		coIt.Set(0.0);

		++mstIt;
		++coIt;

		continue;
	}

	mstExtract->SetExtractionRegion(curMstRegion);
	slvExtract->SetExtractionRegion(curSlvRegion);

    ImageType::SizeType paddsize;
    paddsize[0] = patchSizeRange / 2;
	paddsize[1] = patchSizeAzimuth / 2;

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
  
    //MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
    //minMax->SetImage(modulus->GetOutput());
    //minMax->ComputeMaximum();
/*
	mean->SetInputImage(modulus->GetOutput());
	IndexType windowIndex;
	windowIndex[0] = patchSizeRange / 2;
	windowIndex[1] = patchSizeAzimuth / 2;
    coIt.Set(mean->EvaluateAtIndex(windowIndex));
*/
	coIt.Set(modulus->GetOutput()->GetPixel(mstIndex));

    std::cout << "Coherence: " << coIt.Get() << std::endl;
  }

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

  //otb::StandardWriterWatcher watcher1(writer,extract,"Extracting from master image.");

  writer->Update();

  return EXIT_SUCCESS;
}

