#ifndef __otbSARCoregistrationImageFilter_h
#define __otbSARCoregistrationImageFilter_h

#include "itkImageToImageFilter.h"
#include "otbComplexInterpolateImageFunction.h"
#include "itkContinuousIndex.h"

#include "otbGenericRSTransform.h"

#include "itkPointSet.h"
#include "otbGridIntersectionPointSetSource.h"
#include "itkFFTComplexToComplexImageFilter.h"

#include "itkFFTShiftImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"

namespace otb
{

/** \class SARCoregistrationImageFilter
 * \brief Computes local offsets between master and slave images by applying GenericRSTransform.
 *
 * This filter tries to find at each location of the fixed image the corresponding best matching
 * patch of radius set by SetRadius() method in the moving image within a search windows whose radius
 * is defined by SetSearchRadius() method.
 *
 * Once a coarse (pixel wise) offset has been found, this match is further refined using dichotomic search
 * until sub-pixel accuracy given by the SetSubPixelAccuracy() is reached.
 *
 * The filter proposes two outputs: GetOutput() return the image of the correlation maximum at each location, and
 * the GetOutputDeformationField() method returns the corresponding offset.
 *
 * If the UseSpacingOn() flag is used, the output deformation field takes the input image spacing into account.
 * otherwise, the deformation field is expressed in pixels (default is On).
 *
 * This filter accepts master and slave images with different sizes and spacing. Metric and search windows radius
 * are expressed in terms of number of pixels in the master image.
 *
 * It is possible to generate an output correlation map and deformation field at a coarser resolution by setting
 * grid step to value higher than 1 (grid step is expressed in terms of number of fixed image pixels).
 * Default value is 1.
 *
 *
 * \example pending/SARCoregistrationImageFilterExample.cxx
 *
 * \sa      
 * \ingroup IntensityImageFilters, Streamed
 */
template <class TInputImage, class TInterpolateFunction, class T0utputCorrelation, class TOutputDeformationField>
class ITK_EXPORT SARCoregistrationImageFilter : public itk::ImageToImageFilter<TInputImage, T0utputCorrelation>
{
public:
  /** Standard class typedefs. */
  typedef SARCoregistrationImageFilter                             Self;
  typedef itk::ImageToImageFilter<TInputImage, T0utputCorrelation> Superclass;
  typedef itk::SmartPointer<Self>                                 Pointer;
  typedef itk::SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SARCoregistrationImageFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef typename T0utputCorrelation::RegionType                 OutputImageRegionType;
  typedef typename TOutputDeformationField::PixelType             DeformationValueType;
  typedef typename TInputImage::Pointer                           InputImagePointerType;
  typedef typename TInputImage::RegionType                        InputImageRegionType;
  typedef typename TInputImage::SizeType                          SizeType;
  typedef typename TInputImage::IndexType                         IndexType;
  typedef typename TInputImage::SpacingType                       SpacingType;
  typedef typename TInputImage::PointType                         PointType;
  typedef typename TInputImage::OffsetType                        OffsetType;
  typedef typename TInputImage::PixelType						  PixelType;
  typedef itk::ConstantBoundaryCondition< TInputImage >           BoundaryConditionType;
  typedef otb::ComplexInterpolateImageFunction<TInputImage, 
		TInterpolateFunction, BoundaryConditionType, double>      InterpolatorType;
  typedef typename InterpolatorType::Pointer                      InterpolatorPointerType;
  typedef itk::ContinuousIndex<double, 2>                         ContinuousIndexType;
  typedef typename itk::Transform<double, 2, 2>                   TransformType;
  typedef typename TransformType::Pointer                         TransformPointerType;

  typedef itk::PointSet< PixelType, TInputImage::ImageDimension > PointSetType;
  typedef otb::GridIntersectionPointSetSource< PointSetType >     PointSetSourceType;
  typedef PointSetType::PointsContainer                           PointsContainerType;
  typedef PointSetType::PointType                                 PointType;

  typedef itk::FFTComplexToComplexImageFilter< 
	  PixelType::value_type, ImageType::ImageDimension >		  FFTType;
  typedef FFTType::OutputImageType                                FFTOutputImageType;
  typedef itk::FFTShiftImageFilter<FFTOutputImageType,FFTOutputImageType> ShiftFilterType;
  typedef itk::ConstantPadImageFilter<TInputImage, TInputImage>   PadFilterType;
  typedef otb::Image<double,2>                                    RealImageType;
  typedef itk::ComplexToModulusImageFilter<FFTOutputImageType,
	  RealImageType>											  ModulusFilterType;
  typedef itk::DivideImageFilter<FFTOutputImageType,
	  RealImageType,FFTOutputImageType>							  DivideFilterType;
  typedef itk::MinimumMaximumImageCalculator<RealImageType>       MinMaxCalculatorType;

  /** Set/Get the interpolator used to interpolate moving image at non-grid positions */
  itkSetObjectMacro(Interpolator, InterpolatorType);
  itkGetObjectMacro(Interpolator, InterpolatorType);

  /** Connect one of the operands for registration */
  void SetMasterInput( const TInputImage * image);

  /** Connect one of the operands for registration */
  void SetSlaveInput( const TInputImage * image);

  /** Get the inputs */
  const TInputImage * GetMasterInput();
  const TInputImage * GetSlaveInput();

  /** Get the output deformation field */
  TOutputDeformationField * GetOutputDeformationField();
 
  /** Set the number of desired initial tie points in each dimension */
  itkSetMacro(TiePointsPerDim, unsigned int);
  itkGetMacro(TiePointsPerDim, unsigned int);

  /** Set the size of the area on which correlation is computed */
  itkSetMacro(PatchSizePerDim, SizeType);
  itkGetMacro(PatchSizePerDim, SizeType);

  /** Set the searh radius for fine registration */
  itkSetMacro(SearchRadius, SizeType);
  itkGetMacro(SearchRadius, SizeType);
  
  /** Set/Get subpixel accuracy */
  itkSetMacro(SubPixelAccuracy, double);
  itkGetMacro(SubPixelAccuracy, double);

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

  /** Set/Get the transform for the initial offset */
  itkSetObjectMacro(Transform, TransformType);
  itkGetConstObjectMacro(Transform, TransformType);

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
  PixelType						m_CorrelationThreshold;

  /** The number of tie points */
  unsigned int					m_TiePointsPerDim;

  /** The size for correlation */
  SizeType                      m_PatchSizePerDim;

  /** The search radius */
  SizeType                      m_SearchRadius;

  /** If true, deformation field uses spacing. Otherwise, uses pixel grid */
  bool                          m_UseSpacing;

  /** Search step */
  double                        m_SubPixelAccuracy;

  /** The interpolator */
  InterpolatorPointerType       m_Interpolator;

  /** The translation */
  TranslationPointerType        m_Translation;

  /** Grid step */
  OffsetType                    m_GridStep;

  /** Transform for initial offset */
  TransformPointerType          m_Transform;

};

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbSARCoregistrationImageFilter.txx"
#endif

#endif
