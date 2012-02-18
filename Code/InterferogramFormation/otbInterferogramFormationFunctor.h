#ifndef __otbInterferogramFormationFunctor_h
#define __otbInterferogramFormationFunctor_h

namespace otb {

namespace Functor {

/**
 * \class InterferogramFormationFunctor
 * \brief Function object which computes interferogram from two coregistered images.
 *
 */
template<class TInput1, class TInput2, class TOutput>
class InterferogramFormationFunctor
{
public:
  InterferogramFormationFunctor() {};
  ~InterferogramFormationFunctor() {}

  inline TOutput operator()( const TInput1 & itA, const TInput2 & itB)
    {
      TOutput result = itk::NumericTraits<TOutput>::Zero;
      double normA = 0.0;
      double normB = 0.0;
      for (unsigned long pos = 0; pos< itA.Size(); ++pos)
      {
        //TODO: add multiplication by G(m,i) function
        result += static_cast<TOutput>(itA.GetPixel(pos)* std::conj(itB.GetPixel(pos)));
        normA += itA.GetPixel(pos).real()*itA.GetPixel(pos).real()
            + itA.GetPixel(pos).imag()*itA.GetPixel(pos).imag();
        normB += itB.GetPixel(pos).real()*itB.GetPixel(pos).real()
            + itB.GetPixel(pos).imag()*itB.GetPixel(pos).imag();

      }
      if ((normA != 0) && (normB != 0))
        return static_cast<TOutput>( result/ (vcl_sqrt(normA)*vcl_sqrt(normB)) );
      else
        return itk::NumericTraits<TOutput>::Zero;
    }
};

} // end namespace functor

} // end namespace otb

#endif
