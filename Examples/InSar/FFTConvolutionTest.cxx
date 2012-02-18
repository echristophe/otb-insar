#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbStreamingImageFileWriter.h"
#include "otbImageFileWriter.h"

#include "itkFFTComplexToComplexImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"

typedef std::complex< double >                                         PixelType;
typedef otb::Image< PixelType,2 >                                      ImageType;
typedef otb::Image<double,2>                                           RealImageType;
typedef itk::FFTComplexToComplexImageFilter< PixelType::value_type, 
                                           ImageType::ImageDimension > FFTType;
typedef FFTType::OutputImageType                                       FFTOutputImageType;
typedef itk::FFTShiftImageFilter<FFTOutputImageType,FFTOutputImageType> ShiftFilterType;
typedef itk::ConstantPadImageFilter<ImageType,ImageType>          PadFilterType;
typedef FFTType::TransformDirectionType                                FFTDirectionType;
typedef otb::ImageFileReader< ImageType >                              ReaderType;
typedef otb::ImageFileWriter<FFTOutputImageType>                       FFTWriterType;
typedef otb::StreamingImageFileWriter< ImageType >                     WriterType
;typedef otb::StreamingImageFileWriter< RealImageType >                RealWriterType;
typedef itk::ComplexToModulusImageFilter<FFTOutputImageType,RealImageType>      ModulusFilterType;
typedef itk::MinimumMaximumImageCalculator<RealImageType>              MinMaxCalculatorType;
typedef itk::DivideImageFilter<FFTOutputImageType,RealImageType,FFTOutputImageType> DivideFilterType;
class ComplexConjugateProduct
{
public:
  inline PixelType operator()(const PixelType & a, const PixelType & b) const
  {
    return a * vcl_conj(b);
  }
};

typedef itk::BinaryFunctorImageFilter<FFTOutputImageType,FFTOutputImageType,FFTOutputImageType,ComplexConjugateProduct> ConjugateProductFilterType;

int main(int argc, char * argv[])
{
  if(argc!=3)
    {
    std::cerr<<argv[0]<<" infname1 infname2"<<std::endl;
    return EXIT_FAILURE;
    }

  const char * infname1 = argv[1];
  const char * infname2 = argv[2];

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName(infname1);
  reader1->UpdateOutputInformation();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName(infname2);
  reader2->UpdateOutputInformation();

  if(reader1->GetOutput()->GetLargestPossibleRegion() 
     != reader2->GetOutput()->GetLargestPossibleRegion())
    {
    std::cerr<<"Images sizes must match."<<std::endl;
    return EXIT_FAILURE;
    }
  
  ImageType::SizeType size = reader1->GetOutput()->GetLargestPossibleRegion().GetSize();
  size[0]=size[0]/2;
  size[1]=size[1]/2;


  PadFilterType::Pointer pad1 = PadFilterType::New();
  pad1->SetInput(reader1->GetOutput());
  pad1->SetPadBound(size);
  std::cout<<"Pad size: "<<size<<std::endl;

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetFileName("paded1.tif");
  writer1->SetInput(pad1->GetOutput());
  writer1->Update();

  PadFilterType::Pointer pad2 = PadFilterType::New();
  pad2->SetInput(reader2->GetOutput());
  pad2->SetPadBound(size);

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName("paded2.tif");
  writer2->SetInput(pad2->GetOutput());
  writer2->Update();

  FFTType::Pointer fft = FFTType::New();
  fft->SetInput(pad1->GetOutput());
  fft->Update();

  FFTOutputImageType::Pointer fft1 = fft->GetOutput();

  FFTWriterType::Pointer writer3 = FFTWriterType::New();
  writer3->SetFileName("fft1.tif");
  writer3->SetInput(fft1);
  writer3->Update();

  fft = FFTType::New();
  fft->SetInput(pad2->GetOutput());
  fft->Update();

  FFTOutputImageType::Pointer fft2 = fft->GetOutput();


  FFTWriterType::Pointer writer4 = FFTWriterType::New();
  writer4->SetFileName("fft2.tif");
  writer4->SetInput(fft2);
  writer4->Update();

  ConjugateProductFilterType::Pointer conjProd = ConjugateProductFilterType::New();
  conjProd->SetInput1(fft1);
  conjProd->SetInput2(fft2);

  // Try to normalise to get the normalized coherence coeff
  ModulusFilterType::Pointer conjProdMod = ModulusFilterType::New();
  conjProdMod->SetInput(conjProd->GetOutput());
  
  DivideFilterType::Pointer conjProdNorm = DivideFilterType::New();
  conjProdNorm->SetInput1(conjProd->GetOutput());
  conjProdNorm->SetInput2(conjProdMod->GetOutput());

  FFTWriterType::Pointer writer5 = FFTWriterType::New();
  writer5->SetFileName("fftprod.tif");
  writer5->SetInput(conjProdNorm->GetOutput());
  writer5->Update();

  fft = FFTType::New();
  fft->SetInput(conjProdNorm->GetOutput());
  fft->SetTransformDirection(FFTType::INVERSE);
  fft->Update();
  
  FFTOutputImageType::Pointer ifft = fft->GetOutput();

  ShiftFilterType::Pointer shifter = ShiftFilterType::New();
  shifter->SetInput(ifft);

  FFTWriterType::Pointer writer6 = FFTWriterType::New();
  writer6->SetFileName("invfftprod.tif");
  writer6->SetInput(shifter->GetOutput());
  writer6->Update();

  ModulusFilterType::Pointer modulus = ModulusFilterType::New();
  modulus->SetInput(shifter->GetOutput());
  modulus->Update();
  
  RealWriterType::Pointer writer7 = RealWriterType::New();
  writer7->SetInput(modulus->GetOutput());
  writer7->SetFileName("modinvfftprod.tif");
  writer7->Update();

  MinMaxCalculatorType::Pointer minMax = MinMaxCalculatorType::New();
  minMax->SetImage(modulus->GetOutput());
  minMax->ComputeMaximum();

  std::cout<<"Maximum:  "<<minMax->GetMaximum()<<std::endl;
  std::cout<<"Position: "<<minMax->GetIndexOfMaximum()<<std::endl;

  std::cout<<"Estimated offset: "<<std::endl;
  std::cout<<(int)minMax->GetIndexOfMaximum()[0]-(int)reader1->GetOutput()->GetLargestPossibleRegion().GetSize()[0]/2<<"\t";
  std::cout<<(int)minMax->GetIndexOfMaximum()[1]-(int)reader1->GetOutput()->GetLargestPossibleRegion().GetSize()[1]/2<<std::endl;

  return EXIT_SUCCESS;
}

