#ifndef __otbFlatEarthRemovalFunctorFunctor_h
#define __otbFlatEarthRemovalFunctorFunctor_h

namespace otb {

namespace Functor {

/**
 * \class FlatEarthRemovalFunctorFunctor
 * \brief Function object which computes interferogram from two coregistered images.
 *
 */
template<class TInput, class TOutput>
class FlatEarthRemovalFunctor
{
public:
  FlatEarthRemovalFunctor() : m_RangeFrequency(0.0), m_AzimuthFrequency(0.0), m_Phase(0.0) {};
  ~FlatEarthRemovalFunctor() {}

  typedef itk::Index<2> IndexType;

  inline TOutput operator()( const TInput& inPix, IndexType index )
    {
	//TOutput phase(vcl_cos(otb::CONST_2_PI*m_RangeFrequency*index[0]+otb::CONST_2_PI*m_AzimuthFrequency*index[1]+m_Phase), -(vcl_sin(otb::CONST_2_PI*m_RangeFrequency*index[0]+otb::CONST_2_PI*m_AzimuthFrequency*index[1]+m_Phase)));
	//TOutput phase(vcl_cos(otb::CONST_2_PI*0.98*m_RangeFrequency*index[0]), -(vcl_sin(otb::CONST_2_PI*0.98*m_RangeFrequency*index[0]))); // 0.98 is a magic number
	TOutput phase(vcl_cos(otb::CONST_2_PI*m_RangeFrequency*index[0]), -(vcl_sin(otb::CONST_2_PI*m_RangeFrequency*index[0])));
      return static_cast<TOutput>( inPix * phase );
    }
  void SetRangeFrequency(double value)
  {
    m_RangeFrequency = value;
  }

  void SetAzimuthFrequency(double value)
  {
    m_AzimuthFrequency = value;
  }

  void SetPhaseAtMaxFreq(double value)
  {
	  m_Phase = value;
  }

  private:
	  double m_RangeFrequency;
	  double m_AzimuthFrequency;
	  double m_Phase;
};

} // end namespace functor

} // end namespace otb

#endif
