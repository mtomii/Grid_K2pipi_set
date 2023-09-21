/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AFourQuarkContractionMT.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AFourQuarkContractionMT_hpp_
#define Hadrons_MContraction_A2AFourQuarkContractionMT_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  From closed loop from all-to-all vectors                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AFourQuarkContractionMTPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFourQuarkContractionMTPar,
				    int, nt,
				    std::string, sctypes,
                                    std::string, output,
                                    std::string, mat1,
				    std::string, mat2,
				    std::string, gammas1,
				    std::string, gammas2);
};

template <typename FImpl>
class TA2AFourQuarkContractionMT: public Module<A2AFourQuarkContractionMTPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2AFourQuarkContractionMT(const std::string name);
    // destructor
    virtual ~TA2AFourQuarkContractionMT(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::vector<std::vector<Gamma::Algebra> >       gamma1_;
  std::vector<std::vector<Gamma::Algebra> >       gamma2_;
  std::vector<std::string> nameg1_;
  std::vector<std::string> nameg2_;
  std::vector<int> types_;
  std::string mat1_;
  std::string mat2_;
};

MODULE_REGISTER_TMP(A2AFourQuarkContractionMT, TA2AFourQuarkContractionMT<FIMPL>, MContraction);

/******************************************************************************
 *                 TA2AFourQuarkContractionMT implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AFourQuarkContractionMT<FImpl>::TA2AFourQuarkContractionMT(const std::string name)
: Module<A2AFourQuarkContractionMTPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionMT<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().mat1, par().mat2};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionMT<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContractionMT<FImpl>::setup(void)
{
  gamma1_.clear();
  gamma2_.clear();
  std::vector<std::string> tmp1 = strToVec<std::string>(par().gammas1);
  std::vector<std::string> tmp2 = strToVec<std::string>(par().gammas2);
  nameg1_ = tmp1;
  nameg2_ = tmp2;
  assert( tmp1.size() == tmp2.size() );
  for ( int ig = 0; ig < tmp1.size(); ++ig ) {
    std::vector<Gamma::Algebra> vec;
    vec.clear();
    if ( tmp1[ig] == "GammaMU" ) {
      vec = {
	Gamma::Algebra::GammaX,
	Gamma::Algebra::GammaY,
	Gamma::Algebra::GammaZ,
	Gamma::Algebra::GammaT
      };
    } else if ( tmp1[ig] == "GammaMUGamma5" ) {
      vec = {
	Gamma::Algebra::GammaXGamma5,
	Gamma::Algebra::GammaYGamma5,
	Gamma::Algebra::GammaZGamma5,
	Gamma::Algebra::GammaTGamma5
      };
    } else {
      vec = strToVec<Gamma::Algebra>(tmp1[ig]);
    }
    gamma1_.push_back(vec);

    vec.clear();
    if ( tmp2[ig] == "GammaMU" ) {
      vec = {
	Gamma::Algebra::GammaX,
	Gamma::Algebra::GammaY,
	Gamma::Algebra::GammaZ,
	Gamma::Algebra::GammaT
      };
    } else if ( tmp2[ig] == "GammaMUGamma5" ) {
      vec = {
	Gamma::Algebra::GammaXGamma5,
	Gamma::Algebra::GammaYGamma5,
	Gamma::Algebra::GammaZGamma5,
	Gamma::Algebra::GammaTGamma5
      };
    } else {
      vec = strToVec<Gamma::Algebra>(tmp2[ig]);
    }
    gamma2_.push_back(vec);
  }
  types_ = strToVec<int>(par().sctypes);
  mat1_ = par().mat1;
  mat2_ = par().mat2;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  //envCreateLat(PropagatorField, getName());
  envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContractionMT<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Scalar_v;

  auto &mat1 = envGet(std::vector<SpinColourMatrix_v>, par().mat1);
  auto &mat2 = envGet(std::vector<SpinColourMatrix_v>, par().mat2);
  auto &nt   = envGet(int, par().nt);

  int vol3d = mat1.size() / nt;
  assert ( mat1.size() == mat2.size() );

  std::vector<Scalar_v> corr0(nt,Zero());
  std::vector<std::vector<Scalar_v> > corr(nt,Zero());
  int num_corr = gamma1_.size() * types_.size();
  corr.assign(num_corr,corr0);

  /*  
  thread_for(ix,mat1.size(),{
    for(int ix=0;ix<MFrvol;ix++){
      int sv = ix+vol3d*it;
      for(int s1=0;s1<Ns;s1++)
      for(int s2=0;s2<Ns;s2++)
      for(int c1=0;c1<Nc;c1++)
      for(int c2=0;c2<Nc;c2++){
	mat[ix]()(s1,s2)(c1,c2) += leftV[sv]()(s1)(c1) * rightW[sv]()(s2)(c2);
      }
    }
  });
  */
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AFourQuarkContractionMT_hpp_
