#ifndef __otbSARCoregistrationImageFilter_h
#define __otbSARCoregistrationImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkContinuousIndex.h"
#include "otbExtractROI.h"

#include "otbGenericRSTransform.h"
#include "otbLeastSquareAffineTransformEstimator.h"

#include "otbComplexInterpolateImageFunction.h"
#include "otbImageNormalizeZeroFrequencyCalculator.h"

#include "otbStreamingResampleImageFilter.h"

#include "itkPointSet.h"
#include "otbGridIntersectionPointSetSource.h"
#include "itkFFTComplexToComplexImageFilter.h"

#include "itkFFTShiftImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"

//#include "itkLinearInterpolateImageFunction.h"
//#include "otbBSplinesInterpolateDeformationFieldGenerator.h"

namespace otb
{

/** Complex conjugate product functor class */
class ComplexConjugateProduct
{
public:
	typedef std::complex< double > ValueType;
  inline ValueType operator()(const ValueType & a, const ValueType & b) const
  {
    return a * vcl_conj(b);
  }
};

/** \class SARCoregistrationImageFilter
 * \brief Computes local offsets between master and slave images by applying GenericRSTransform
 * and coarse cross-correlation with affinetransform.
 *
 * This filter tries to find at each grid location of the master image the corresponding best matching
 * patch in the slave image.
 *
 * This filter accepts master and slave images with different sizes and spacing. Metric and search windows radius
 * are expressed in terms of number of pixels in the master image.
 *
 *
 * \example pending/SARCoregistrationImageFilterExample.cxx
 *
 * \sa      
 * \ingroup IntensityImageFilters, Streamed
 */
template <class TInputImage, class TInterpolateFunction>
class ITK_EXPORT SARCoregistrationImageFilter : public itk::ImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard class typedefs. */
  typedef SARCoregistrationImageFilter                          Self;
  typedef itk::ImageToImageFilter<TInputImage, TInputImage>		Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SARCoregistrationImageFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef typename TInputImage::Pointer                           InputImagePointerType;
  typedef typename TInputImage::RegionType                        InputImageRegionType;
  typedef typename TInputImage::SizeType                          SizeType;
  typedef typename TInputImage::IndexType                         IndexType;
  typedef typename TInputImage::SpacingType                       SpacingType;
  typedef typename TInputImage::RegionType						  RegionType;
  typedef typename TInputImage::OffsetType                        OffsetType;
  typedef typename TInputImage::PixelType						  PixelType;
  typedef itk::ContinuousIndex<double, 
									TInputImage::ImageDimension>  ContinuousIndexType;
  
  typedef typename otb::ExtractROI< PixelType, PixelType >        ExtractFilterType;

  typedef typename itk::PointSet< PixelType, 
									TInputImage::ImageDimension > PointSetType;
  typedef typename otb::GridIntersectionPointSetSource< 
												PointSetType >    PointSetSourceType;
  typedef typename PointSetType::PointsContainer                  PointsContainerType;
  typedef typename PointSetType::PointType                        PointType;

  typedef typename otb::GenericRSTransform<  >					  TransformType;
  typedef typename TransformType::Pointer						  TransformPointerType;
  typedef typename itk::Point< typename PixelType::value_type, 
									TInputImage::ImageDimension > LSQPointType;
  typedef typename otb::LeastSquareAffineTransformEstimator< 
												LSQPointType >    EstimateFilterType;
  typedef typename TInterpolateFunction							  FunctionType;
  typedef typename itk::ConstantBoundaryCondition< TInputImage >  BoundaryConditionType;
  typedef typename PixelType::value_type						  CoordRepType;
  typedef typename otb::ComplexInterpolateImageFunction< 
				TInputImage, FunctionType, BoundaryConditionType, 
												CoordRepType >	  InterpolateType;
  typedef typename InterpolateType::Pointer						  InterpolatePointerType;
  typedef typename otb::ImageNormalizeZeroFrequencyCalculator<
													TInputImage > NormalizeZeroFrequencyType;

  typedef typename itk::FFTComplexToComplexImageFilter< typename
	  PixelType::value_type, TInputImage::ImageDimension >		  FFTType;
  typedef typename FFTType::OutputImageType                       FFTOutputImageType;
  typedef typename itk::FFTShiftImageFilter< FFTOutputImageType,
											FFTOutputImageType >  ShiftFilterType;
  typedef typename itk::ConstantPadImageFilter< TInputImage, 
													TInputImage > PadFilterType;
  typedef typename otb::Image<double,2>                           RealImageType;
  typedef typename itk::ComplexToModulusImageFilter< 
							FFTOutputImageType,	RealImageType >	  ModulusFilterType;
  typedef typename itk::DivideImageFilter< FFTOutputImageType,
							RealImageType,FFTOutputImageType >	  DivideFilterType;
  typedef typename itk::MinimumMaximumImageCalculator<
												RealImageType >   MinMaxCalculatorType;

  typedef typename otb::StreamingResampleImageFilter< 
									TInputImage, TInputImage >    ResampleFilterType;

  typedef itk::BinaryFunctorImageFilter<FFTOutputImageType,FFTOutputImageType,FFTOutputImageType,ComplexConjugateProduct> ConjugateProductFilterType;

  /** Set/Get the interpolator used to interpolate moving image at non-grid positions */
  itkSetObjectMacro(Interpolator, InterpolateType);
  itkGetObjectMacro(Interpolator, InterpolateType);

  /** Connect one of the operands for registration */
  void SetMasterInput( const TInputImage * image);

  /** Connect one of the operands for registration */
  void SetSlaveInput( const TInputImage * image);

  /** Get the inputs */
  const TInputImage * GetMasterInput();
  const TInputImage * GetSlaveInput();
 
  /** Set the coarse cross-correlation threshold */
  itkSetMacro(CorrelationThreshold, typename PixelType::value_type);
  itkGetMacro(CorrelationThreshold, typename PixelType::value_type);

  /** Set the number of desired initial tie points in each dimension */
  itkSetMacro(TiePointsPerDim, unsigned int);
  itkGetMacro(TiePointsPerDim, unsigned int);

  /** Set the size of the area on which correlation is computed */
  itkSetMacro(PatchSizePerDim, unsigned int);
  itkGetMacro(PatchSizePerDim, unsigned int);

  /** Set the searh radius for fine registration */
  itkSetMacro(SearchRadius, SizeType);
  itkGetMacro(SearchRadius, SizeType);
  
  /** Set/Get subpixel accuracy */
  itkSetMacro(SubPixelAccuracy, double);
  itkGetMacro(SubPixelAccuracy, double);

  /** Set/Get Perform fine registration flag */
  itkSetMacro(PerformFine, bool);
  itkGetConstMacro(PerformFine, bool);

  /** Set/Get fine coherency threshold */
  itkSetMacro(CoherencyThreshold, typename PixelType::value_type);
  itkGetMacro(CoherencyThreshold, typename PixelType::value_type);

  /** Set/Get fine window size */
  itkSetMacro(CoherencyWindowSizePerDim, unsigned int);
  itkGetMacro(CoherencyWindowSizePerDim, unsigned int);

  /** Set/Get dem directory */
  itkSetMacro(DEMDir, char *);
  itkGetMacro(DEMDir, const char *);

  /** Set/Get Use DEM flag */
  itkSetMacro(UseDEM, bool);
  itkGetConstMacro(UseDEM, bool);

  /** True if deformation field takes spacing into account. False otherwise */
  itkSetMacro(UseSpacing, bool);
  itkBooleanMacro(UseSpacing);

  /** Set the grid step */
  itkSetMacro(GridStep, OffsetType);
  itkGetConstReferenceMacro(GridStep, OffsetType);

 /** Set unsigned int radius */
  void SetSearchRadius(unsigned int radius)
  {
    m_SearchRadius.Fill(radius);
  }

  /** Set unsigned int grid step */
  void SetGridStep(unsigned int step)
  {
    m_GridStep.Fill(step);
  }

protected:
  /** Constructor */
  SARCoregistrationImageFilter();
  /** Destructor */
  virtual ~SARCoregistrationImageFilter() {};

  /** Threaded generate data */
  virtual void GenerateData();

  /** Generate the input requested regions  */
  virtual void GenerateInputRequestedRegion(void);

  /** Generate output information */
  virtual void GenerateOutputInformation(void);

private:
  SARCoregistrationImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Coarse correlation threshold */
  typename PixelType::value_type	m_CorrelationThreshold;

  /** The number of tie points */
  unsigned int						m_TiePointsPerDim;

  /** The size for correlation */
  unsigned int						m_PatchSizePerDim;

  /** The search radius */
  SizeType							m_SearchRadius;

  /** If true, deformation field uses spacing. Otherwise, uses pixel grid */
  bool								m_UseSpacing;

  /** Search step */
  double							m_SubPixelAccuracy;

  /** Fine correlation threshold */
  typename PixelType::value_type	m_CoherencyThreshold;

  /** Fine coherency window size */
  unsigned int						m_CoherencyWindowSizePerDim;

  /** Perform fine registration flag */
  bool								m_PerformFine;

  /** The interpolator */
  InterpolatePointerType			m_Interpolator;

  /** Grid step */
  OffsetType						m_GridStep;

  /** Transform for initial offset */
  TransformPointerType				m_Transform;

  /** DEM directory for transform */
  const char *							m_DEMDir;

  bool								m_UseDEM;

};

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbSARCoregistrationImageFilter.txx"
#endif

#endif
