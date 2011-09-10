/*=========================================================================

   Copyright 2011 Patrick IMBO
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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageNormalizeZeroFrequencyCalculator.h"
#include "itkFFTComplexToComplexImageFilter.h"

int 
otbNormalizeZeroFrequencyCalculatorTest(int argc, char * argv[] )
{
	typedef double													InternalPixelType;   
	typedef std::complex<InternalPixelType>							InputPixelType;
	typedef otb::Image<InputPixelType, 2>							ImageType;
	typedef itk::FFTComplexToComplexImageFilter<InternalPixelType,2> FFTImageFilterType;
	typedef otb::ImageNormalizeZeroFrequencyCalculator<FFTImageFilterType::OutputImageType>	NormalizeZeroFreqCalculatorType;

	unsigned int imageSizeX = std::atoi(argv[1]);
	unsigned int imageSizeY = std::atoi(argv[2]);
	float diracX = std::atof(argv[3]);
	float diracY = std::atof(argv[4]);

    /* Define the image size and physical coordinates */
    itk::Size<2> size = {{imageSizeX, imageSizeY}};
    double origin [2] = { 0.0,   0.0};
    double spacing[2] = { 1.0/imageSizeX,   1.0/imageSizeY};

    /* Allocate a simple test image */
    ImageType::Pointer image = ImageType::New();

    ImageType::RegionType region;
    region.SetSize(size);
    image->SetLargestPossibleRegion(region);
    image->SetBufferedRegion(region);
    image->SetRequestedRegion(region);

    /* Set origin and spacing of physical coordinates */
    image->SetOrigin(origin);
    image->SetSpacing(spacing);
    image->Allocate();
    
    image->FillBuffer( itk::NumericTraits<InputPixelType>::Zero );


    /* Define positions of Dirac in index coordinates */
	itk::Index<2> diracIndex = {{diracX*imageSizeX, diracY*imageSizeY}};
	InputPixelType diracValue(1.0,0.0);
	std::cout << "diracIndex : " << diracIndex[0] << "and " << diracIndex[1] <<std::endl;

    image->SetPixel(diracIndex, diracValue);


	/* Generate image*/
	FFTImageFilterType::Pointer inverseFFTImage = FFTImageFilterType::New();

	inverseFFTImage->SetInput(image);
	inverseFFTImage->SetTransformDirection(FFTImageFilterType::INVERSE);
	inverseFFTImage->Update();

	/* Compute Normalize Zero frequency */	
	NormalizeZeroFreqCalculatorType::Pointer zeroFreqCalculator = NormalizeZeroFreqCalculatorType::New();
	zeroFreqCalculator->SetImage(inverseFFTImage->GetOutput());
	zeroFreqCalculator->Compute();

	/* Printout info */
	zeroFreqCalculator->Print( std::cout );

	if( (std::abs(zeroFreqCalculator->GetNormalizeZeroFrequency()[0]-diracX) > spacing[0]) &&
		(std::abs(zeroFreqCalculator->GetNormalizeZeroFrequency()[1]-diracY) > spacing[1])     )
	{
		return EXIT_FAILURE;
	}

    return EXIT_SUCCESS;
}
