#ifndef __otbFlatEarthRemovalImageFilter_h
#define __otbFlatEarthRemovalImageFilter_h

#include "itkImageToImageFilter.h"
#include "otbExtractROI.h"

#include "itkImageRegionIterator.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "otbSubsampledImageRegionConstIterator.h"

#include "otbWindowedSincInterpolateImageBlackmanFunction.h"
#include "otbComplexInterpolateImageFunction.h"
#include "otbImageNormalizeZeroFrequencyCalculator.h"
#include "itkScaleTransform.h"
#include "otbStreamingResampleImageFilter.h"

#include "itkFFTComplexToComplexImageFilter.h"

#include "itkFFTShiftImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "itkMultiplyImageFilter.h"

#include "otbUnaryFunctorWithIndexImageFilter.h"

namespace otb
{

/** Flat Earth Calculation functor class */
template<class TInput, class TOutput>
class FlatEarthPhaseCalculation
{
public:
	typedef itk::Index<2> IndexType;
	inline TOutput operator()(const TInput& inValue, IndexType index)
	{
		//TOutput real = vcl_cos(otb::CONST_2_PI*m_RangeRate*index[0])*vcl_cos(otb::CONST_2_PI*m_AzimuthRate*index[1])*vcl_cos(m_FringePhase);
		//TOutput imag = -(vcl_sin(otb::CONST_2_PI*m_RangeRate*index[0])*vcl_sin(otb::CONST_2_PI*m_AzimuthRate*index[1])*vcl_sin(m_FringePhase);

		TOutput phase(vcl_cos(otb::CONST_2_PI*m_RangeRate*index[0]),
	               -(vcl_sin(otb::CONST_2_PI*m_RangeRate*index[0])));

		return static_cast<TOutput> (phase);
	}

	void SetFringePhase(double value)
	{
		m_FringePhase = value;
	}
	void SetRangeRate(double value)
	{
		m_RangeRate = value;
	}
	void SetAzimuthRate(double value)
	{
		m_AzimuthRate = value;
	}

private:
	double m_FringePhase;
	double m_RangeRate;
	double m_AzimuthRate;
};

/** \class FlatEarthRemovalImageFilter
 * \brief 
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT FlatEarthRemovalImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FlatEarthRemovalImageFilter                          Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage>	Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FlatEarthRemovalImageFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef typename TInputImage::Pointer                           InputImagePointerType;
  typedef typename TInputImage::RegionType                        InputImageRegionType;
  typedef typename TInputImage::SizeType                          SizeType;
  typedef typename TInputImage::IndexType                         IndexType;
  typedef typename TInputImage::SpacingType                       SpacingType;
  typedef typename TInputImage::RegionType						  RegionType;
  typedef typename TInputImage::OffsetType                        OffsetType;
  typedef typename TInputImage::PixelType						  PixelType;
  
  typedef typename otb::ExtractROI< PixelType, PixelType >        ExtractFilterType;

  typedef typename otb::Function::BlackmanWindowFunction< typename PixelType::value_type >	  FunctionType;
  typedef typename itk::ConstantBoundaryCondition< TInputImage >  BoundaryConditionType;
  typedef typename PixelType::value_type						  CoordRepType;
  typedef typename otb::ComplexInterpolateImageFunction< 
				TInputImage, FunctionType, BoundaryConditionType, 
												CoordRepType >	  InterpolateType;
  typedef typename otb::ImageNormalizeZeroFrequencyCalculator< 
													TInputImage > NormalizeFrequencyType;
  typedef typename itk::ScaleTransform<
	typename PixelType::value_type, TInputImage::ImageDimension > TransformType;

  typedef typename otb::StreamingResampleImageFilter< 
									TInputImage, TInputImage >    ResampleFilterType;

  typedef typename itk::FFTComplexToComplexImageFilter< typename
	  PixelType::value_type, TInputImage::ImageDimension >		  FFTType;
  typedef typename FFTType::OutputImageType                       FFTOutputImageType;
  typedef typename itk::ConstantPadImageFilter< TInputImage, 
													TInputImage > PadFilterType;
    typedef typename itk::ComplexToModulusImageFilter< 
							FFTOutputImageType,	TOutputImage >	  ModulusFilterType;
  typedef typename itk::ComplexToPhaseImageFilter< FFTOutputImageType, TOutputImage > PhaseFilterType;
  typedef typename itk::MultiplyImageFilter< TInputImage, TInputImage, TInputImage > MultiplyFilterType;
  typedef typename itk::MinimumMaximumImageCalculator<
												TOutputImage >   MinMaxCalculatorType;

  //typedef typename itk::ImageRegionConstIterator< TInputImage >		  BlockIteratorType;
  typedef typename itk::ImageLinearConstIteratorWithIndex< TInputImage > LinearIteratorType;
  typedef typename otb::SubsampledImageRegionConstIterator< TInputImage> SubsampledRegionIteratorType;
  typedef typename itk::ImageRegionIterator< TOutputImage >				OutRegionIteratorType;
  typedef typename itk::ImageRegionIterator< TInputImage >				InRegionIteratorType;

  typedef otb::FlatEarthPhaseCalculation< typename TOutputImage::PixelType, typename TInputImage::PixelType > FlatEarthPhaseFunctorType;
  typedef otb::UnaryFunctorWithIndexImageFilter< TOutputImage, TInputImage, FlatEarthPhaseFunctorType > FlatEarthPhaseCalculationType;

  /** Set the size of the area on which correlation is computed */
  itkSetMacro(PatchSizePerDim, unsigned int);
  itkGetMacro(PatchSizePerDim, unsigned int);
  
  /** Set/Get interpolation pad size */
  itkSetMacro(PadSizePerDim, double);
  itkGetMacro(PadSizePerDim, double);

protected:
  /** Constructor */
  FlatEarthRemovalImageFilter();
  /** Destructor */
  virtual ~FlatEarthRemovalImageFilter() {};

  /** Threaded generate data */
  virtual void GenerateData();

  /** Generate the input requested regions  */
  virtual void GenerateInputRequestedRegion(void);

  /** Generate output information */
  virtual void GenerateOutputInformation(void);

private:
  FlatEarthRemovalImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The size for correlation */
  unsigned int						m_PatchSizePerDim;

  /** Interpolation size */
  unsigned int						m_PadSizePerDim;
};

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbFlatEarthRemovalImageFilter.txx"
#endif

#endif
