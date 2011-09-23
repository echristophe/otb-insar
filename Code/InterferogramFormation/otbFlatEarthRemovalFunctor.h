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
  FlatEarthRemovalFunctor() : m_RangeFrequency(0.0) {};
  ~FlatEarthRemovalFunctor() {}

  typedef itk::Index<2> IndexType;

  inline TOutput operator()( const TInput& inPix, IndexType index )
    {
      TOutput phase(cos(otb::CONST_2PI*m_RangeFrequency*index[1]), sin(otb::CONST_2PI*m_RangeFrequency*index[1]));
      return static_cast<TOutput>( inPix * phase );
    }
  void SetRangeFrequency(double value)
  {
    m_RangeFrequency = value;
  }

  private:
	  double m_RangeFrequency;
};

} // end namespace functor

} // end namespace otb

#endif
