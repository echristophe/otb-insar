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
#ifndef __otbMultivariateRationalTransform_h
#define __otbMultivariateRationalTransform_h

#include "otbTransform.h"
#include "itkMacro.h"
#include "itkIndex.h"
#include "itkVectorContainer.h"

namespace otb
{
/** \class MultivariateRationalTransform
 *  \brief This class implements a Multivariate rational transfom
 *
 *  A Multivariate rational transform is a quotient of two polynomial functions.
 *
 *  The degree of the numerator and denominator polynomial functions
 *  can be set using the appropriate setters.
 *
 *  The number of parameters is then the number of dimensions times
 *  the numerator degree plus one times the denominator degree plus
 *  one.
 *
 *  Parameters in the parameters vector are in the following order:
 *  dim0num0 dim0num1 ... dim0numN dim0denom0 dim0denom1
 *  ... dim0denomM ... dim1num0 ... dimDdenomM.
 *
 * \ingroup Transform
 **/

template <class TScalarType = double,
          unsigned int Dimension = 2>
class ITK_EXPORT MultivariateRationalTransform : public Transform<TScalarType, Dimension, Dimension>
{
public:
  /** Standard class typedefs. */
  typedef Transform<TScalarType, Dimension,
    Dimension>                                      Superclass;
  typedef MultivariateRationalTransform             Self;
  typedef itk::SmartPointer<Self>                   Pointer;
  typedef itk::SmartPointer<const Self>             ConstPointer;

  typedef typename Superclass::ScalarType           ScalarType;
  typedef itk::Point<ScalarType, Dimension>         InputPointType;
  typedef itk::Point<ScalarType, Dimension>         OutputPointType;

  typedef itk::Vector<double, Dimension>      SpacingType;
  typedef itk::Point<double, Dimension>       OriginType;
  typedef itk::Index<Dimension>               IndexType;
  typedef itk::VectorContainer<unsigned int, IndexType>   IndexContainerType;
  typedef typename IndexContainerType::Pointer            IndexContainerPointer;

  typedef typename Superclass::InverseTransformBasePointer InverseTransformBasePointer;
  typedef typename Superclass::ParametersType              ParametersType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultivariateRationalTransform, Transform);

  itkStaticConstMacro(SpaceDimension, unsigned int, Dimension);

  /** Set the numerator degree */
  void SetNumeratorDegree(IndexType value)
  {
    this->m_NumeratorDegree = value;
    this->InitalizeParameters();
  }

  void SetNumeratorDegree(unsigned int value)
  {
	 for(unsigned int idx = 0 ; idx < m_NumeratorDegree->Size(); ++idx)
	 {
		IndexType currentIndex = m_NumeratorDegree->GetElement(idx);
		currentIndex.Fill(value);
		m_NumeratorDegree->SetElement(idx,currentIndex);
	 }
    this->InitalizeParameters();
  }


  /** Get the numerator degree */
  //itkGetObjectMacro(NumeratorDegree, IndexContainerType);

  /** Set the numerator degree */
  void SetDenominatorDegree(IndexType value)
  {
    this->m_DenominatorDegree = value;
    this->InitalizeParameters();
  }

  void SetDenominatorDegree(unsigned int value)
  {
	 for(unsigned int idx = 0 ; idx < m_DenominatorDegree->Size(); ++idx)
	 {
		IndexType currentIndex = m_DenominatorDegree->GetElement(idx);
		currentIndex.Fill(value);
		m_DenominatorDegree->SetElement(idx,currentIndex);
	 }
    this->InitalizeParameters();
  }

  /** Get the denominator degree */
//  itkGetConstMacro(DenominatorDegree, IndexType);

  /** The transform point method */
  virtual OutputPointType TransformPoint(const InputPointType& point) const
  {
    // Check for consistency
    if(this->GetNumberOfParameters() != this->m_Parameters.size())
      {
      {
      itkExceptionMacro(<<"Wrong number of parameters: found "<<this->m_Parameters.Size()<<", expected "<<this->GetNumberOfParameters());
      }
      }

    InputPointType  inputPoint = point;
    OutputPointType outputPoint;

    unsigned int coeffOffset = 0;
      
	for(unsigned int dim = 0; dim < SpaceDimension; ++dim)
      {
		  std::vector<double> coefficient(1);
		  coefficient[0]=1.0;
		  for(unsigned int idx = SpaceDimension ; idx > 0  ; idx--)
		  {
			std::vector<double> coeffTemp;
			coeffTemp = coefficient;
			coefficient.clear();
			for(unsigned int deg = 0 ; deg <=m_NumeratorDegree->GetElement(dim)[idx-1] ; deg++)
			{						
				coefficient.insert(coefficient.end(),coeffTemp.begin(),coeffTemp.end());
				for(unsigned int idx2 =0 ; idx2 < coeffTemp.size(); idx2++)
				{
					coeffTemp[idx2] = inputPoint[dim] * coeffTemp[idx2];			
				}
			}
		  }

		  double numerator = 0.0;
		  for(unsigned int idx = 0 ; idx < coefficient.size()  ; ++idx)
		  {
			  numerator += ( coefficient[idx] * this->m_Parameters[idx + coeffOffset]);
		  }

		  coeffOffset += coefficient.size();

		  coefficient.resize(1);
		  coefficient[0]=1.0;
		  for(unsigned int idx = SpaceDimension ; idx > 0  ; idx--)
		  {
			std::vector<double> coeffTemp;
			coeffTemp = coefficient;
			coefficient.clear();
			for(unsigned int deg = 0 ; deg <=m_DenominatorDegree->GetElement(dim)[idx-1] ; ++deg)
			{						
				coefficient.insert(coefficient.end(),coeffTemp.begin(),coeffTemp.end());
				for(unsigned int idx2 =0 ; idx2 < coeffTemp.size(); ++idx2)
				{
					coeffTemp[idx2] = inputPoint[dim] * coeffTemp[idx2];			
				}
			}
		  }

		  double denominator = 0.0;
		  for(unsigned int idx = 0 ; idx < coefficient.size()  ; ++idx)
		  {
			  denominator += ( coefficient[idx] * this->m_Parameters[idx + coeffOffset] );
		  }

		  coeffOffset += coefficient.size();
		  outputPoint[dim] = numerator /denominator;
	}

    // Return the output point
    return outputPoint;
  }
 



  /** Get the number of parameters */
  virtual unsigned int GetNumberOfParameters() const
  {
	  unsigned int numberOfParameters = 0;
      for(unsigned int dim = 0; dim < SpaceDimension; ++dim)
      {
       unsigned int numeratorDegree = 1;
       unsigned int denominatorDegree = 1;
		  for(unsigned int idx = 0 ; idx < SpaceDimension ; ++idx)
		  {
			  numeratorDegree *=(m_NumeratorDegree->GetElement(dim)[idx]+1);
			  denominatorDegree *=(m_DenominatorDegree->GetElement(dim)[idx]+1);
		  }
	  numberOfParameters += (numeratorDegree+denominatorDegree);
	}
	return numberOfParameters;
  }

  // Set parameter method
  virtual void SetParameters(const typename Superclass::ParametersType & params)
  {
    // Check for the appropriate size
    if(params.Size() != this->GetNumberOfParameters())
      {
      itkExceptionMacro(<<"Wrong number of parameters: found "<<params.Size()<<", expected "<<this->GetNumberOfParameters());
      }
    
    // Set parameters
    this->m_Parameters = params;
  }

  /** Initialize Parameters size  */
  void InitalizeParameters()
  {
    this->m_Parameters.SetSize(this->GetNumberOfParameters());
    this->m_Parameters.Fill(0);
    unsigned int dimensionStride = 0;
    for(unsigned int dim = 0; dim < SpaceDimension; ++dim)
      {
       unsigned int numeratorDegree = 1;
       unsigned int denominatorDegree = 1;
		  for(unsigned int idx = 0 ; idx < SpaceDimension ; ++idx)
		  {
			  numeratorDegree *=(m_NumeratorDegree->GetElement(dim)[idx]+1);
			  denominatorDegree *=(m_DenominatorDegree->GetElement(dim)[idx]+1);
		  }
      this->m_Parameters[ dimensionStride + numeratorDegree] = 1.;
	  dimensionStride += (numeratorDegree+denominatorDegree);
	}
  }


protected:
  MultivariateRationalTransform() : 
		Superclass(SpaceDimension, SpaceDimension*4)
    {
		m_NumeratorDegree = IndexContainerType::New();
		m_DenominatorDegree = IndexContainerType::New();
		m_NumeratorDegree->Initialize();
		m_DenominatorDegree->Initialize();
		IndexType currentIndex;
		currentIndex.Fill(0);
		for(unsigned int idx = 0 ; idx < SpaceDimension ; ++idx)
		{
			m_NumeratorDegree->InsertElement(idx,currentIndex);
			m_DenominatorDegree->InsertElement(idx,currentIndex);
		}
		this->InitalizeParameters();
    }


  virtual ~MultivariateRationalTransform() {}

  void PrintSelf(std::ostream& os, itk::Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    //os << indent << "Numerator Degree : " << m_NumeratorDegree << std::endl;
    //os << indent << "Denominator Degree : " << m_DenominatorDegree << std::endl;

  }

private:
  MultivariateRationalTransform(const Self &);    //purposely not implemented
  void operator =(const Self&);    //purposely not implemented

  // Degree of numerator
  IndexContainerPointer m_NumeratorDegree;

  // Degree of denominator
  IndexContainerPointer m_DenominatorDegree;
};

} // namespace otb

#endif
