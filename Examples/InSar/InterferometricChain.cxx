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

#include "otbSARCoregistrationImageFilter.h"

#include "otbGenericRSTransform.h"
#include "otbLeastSquareAffineTransformEstimator.h"

#include "otbWindowedSincInterpolateImageBlackmanFunction.h"
#include "otbComplexInterpolateImageFunction.h"

#include "otbStandardWriterWatcher.h"

#include "otbInterferogramFormationFunctor.h"
#include "otbBinaryFunctorNeighborhoodImageFilter.h"

#include "otbFlatEarthRemovalFunctor.h"

int main(int argc, char* argv[])
{

	if (argc != 19)
    {
		std::cerr << "Usage: " << argv[0] << " masterImage slaveImage tiePointsPerDim patchSizePerDim dem startx starty sizex sizey outfname1 outfname2 coherency_threshold subPixelAccuracy subSampleFactor windowSize fine_threshold outfname3 outfname4";
		return EXIT_FAILURE;
    }

	const unsigned int Dimension = 2;

	typedef std::complex< double >             PixelType;
	typedef otb::Image< PixelType,Dimension >  ImageType;

	const char * masterImage = argv[1];
	const char * slaveImage = argv[2];
	
	unsigned int tiePointsPerDim = atoi(argv[3]);
	unsigned int patchSizePerDim = atoi(argv[4]);
	char * demdir  = argv[5];
	
	ImageType::IndexType index;
	index[0] = atoi(argv[6]);
	index[1] = atoi(argv[7]);

	ImageType::SizeType size;
	size[0] = atoi(argv[8]);
	size[1] = atoi(argv[9]);
	
	const char * outfname1 = argv[10];
	const char * outfname2 = argv[11];
	
	double correlation_threshold = atof(argv[12]);
	double subPixelAccuracy = atof(argv[13]);
	ImageType::IndexType::IndexValueType subSampleFactor = atoi(argv[14]);
	unsigned int windowSizePerDim = atoi(argv[15]);
	double coherency_threshold = atof(argv[16]);
	
	const char * outfname3 = argv[17];
	const char * outfname4 = argv[18];

	typedef otb::ImageFileReader< ImageType > ReaderType;
  
	ReaderType::Pointer master = ReaderType::New();
	ReaderType::Pointer slave = ReaderType::New();

	master->SetFileName(masterImage);
	slave->SetFileName(slaveImage);

	master->GenerateOutputInformation();
	slave->GenerateOutputInformation();

	typedef otb::ExtractROI< PixelType, PixelType > ExtractFilterType;

	ExtractFilterType::Pointer mstExtract = ExtractFilterType::New();
	mstExtract->SetInput(master->GetOutput());

	ImageType::RegionType region;
	region.SetIndex(index);
	region.SetSize(size);

	mstExtract->SetExtractionRegion(region);

	ExtractFilterType::Pointer slvExtract = ExtractFilterType::New();
	slvExtract->SetInput(slave->GetOutput());
	slvExtract->SetExtractionRegion(region);
  

	typedef otb::Function::BlackmanWindowFunction< PixelType::value_type > FunctionType;
    typedef otb::SARCoregistrationImageFilter< ImageType, FunctionType >   CoregistrationFilterType;

	CoregistrationFilterType::Pointer coregistration = CoregistrationFilterType::New();
	coregistration->SetMasterInput(mstExtract->GetOutput());
	coregistration->SetSlaveInput(slvExtract->GetOutput());
	coregistration->SetTiePointsPerDim(tiePointsPerDim);
	coregistration->SetPatchSizePerDim(patchSizePerDim);
	coregistration->SetCorrelationThreshold(correlation_threshold);
	coregistration->SetUseDEM(false);
	//coregistration->SetDEMDir(demdir);
	coregistration->SetPerformFine(false);
	//coregistration->SetSubPixelAccuracy(subPixelAccuracy);
	//coregistration->SetCoherencyThreshold(coherency_threshold);
	//coregistration->SetCoherencyWindowSizePerDim(windowSizePerDim);
  
	typedef otb::StreamingImageFileWriter< ImageType > WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outfname1);
	writer->SetInput(mstExtract->GetOutput());
	writer->Update();

	writer = WriterType::New();
	writer->SetFileName(outfname2);
	writer->SetInput(coregistration->GetOutput());

	otb::StandardWriterWatcher watcher1(writer, coregistration, "Generating coregistered slave image.");
	writer->SetNumberOfDivisionsStrippedStreaming(1);
	writer->Update();

	// Compute interferogram
	typedef otb::Functor::InterferogramFormationFunctor< itk::ConstNeighborhoodIterator<ImageType>, itk::ConstNeighborhoodIterator<ImageType>, ImageType::PixelType >  InterferogramCalculatorType;
	typedef otb::BinaryFunctorNeighborhoodImageFilter< ImageType,ImageType,ImageType,InterferogramCalculatorType > InterferogramFilterType;

	InterferogramFilterType::Pointer interferogram = InterferogramFilterType::New();
	interferogram->SetInput1( mstExtract->GetOutput() );
	interferogram->SetInput2( coregistration->GetOutput() );

	ImageType::SizeType radius;
	radius[0] = 1;
	radius[1] = 1;
	interferogram->SetRadius(radius);

	WriterType::Pointer writerInterferogram = WriterType::New();
	writerInterferogram->SetFileName(outfname4);
	writerInterferogram->SetInput(interferogram->GetOutput());
	writerInterferogram->Update();

	// Remove flat earth component
		// TODO

	return EXIT_SUCCESS;
}

