/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/MultMesonField.hpp

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
#ifndef Hadrons_MContraction_MultMesonField_hpp_
#define Hadrons_MContraction_MultMesonField_hpp_

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

class MultMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MultMesonFieldPar,
				    std::string, mf1,
				    std::string, mf2,
                                    std::string, output,
				    int, tsep,
				    int, lsize,
				    int, rsize);
};

class MultMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MultMesonFieldMetadata,
				    std::string, mf1,
				    std::string, mf2);
};

template <typename FImpl>
class TMultMesonField: public Module<MultMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TMultMesonField(const std::string name);
    // destructor
    virtual ~TMultMesonField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::string mf1_;
  std::string mf2_;
};

MODULE_REGISTER_TMP(MultMesonField, TMultMesonField<FIMPL>, MContraction);

/******************************************************************************
 *                     TMultMesonField implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMultMesonField<FImpl>::TMultMesonField(const std::string name)
: Module<MultMesonFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMultMesonField<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().mf1, par().mf2};

  return in;
}

template <typename FImpl>
std::vector<std::string> TMultMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMultMesonField<FImpl>::setup(void)
{
  std::string mf1 = par().mf1;
  std::string mf2 = par().mf2;
  mf1_ = mf1;
  mf2_ = mf2;
  /*
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, 0, Zero());
  */
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMultMesonField<FImpl>::execute(void)
{
  auto &mf1 = envGet(EigenDiskVector<Complex>, par().mf1);
  auto &mf2 = envGet(EigenDiskVector<Complex>, par().mf2);

  int nt         = env().getDim().back();
  int tsep = par().tsep;
  int lsize = par().lsize;
  int rsize = par().rsize;

  HADRONS_A2AM_IO_TYPE init[nt*lsize*rsize];
  A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mat(init,1,1,nt,lsize,rsize);
  for (unsigned int t = 0; t < nt; ++t ) {
    A2AMatrix<Complex> prod;
    const A2AMatrix<Complex> mat1 = mf1[t];
    const A2AMatrix<Complex> mat2 = mf2[(t+tsep)%nt];
    A2AContraction::mul(prod, mat1, mat2);
    for ( unsigned int i = 0; i < lsize; ++i )
    for ( unsigned int j = 0; j < rsize; ++j ) {
      mat(0,0,t,i,j) = prod(i,j);
    }
  }

  std::stringstream ss;
  ss << mf1_ << "_" << tsep << "_" << mf2_;
  std::string filename
    = par().output + "." + std::to_string(vm().getTrajectory()) 
    + "/" + ss.str() + ".h5";
  MultMesonFieldMetadata md;
  md.mf1 = mf1_;
  md.mf2 = mf2_;

  saveA2AMatrixSet(mat,ss.str(),filename,md,lsize,rsize,nt);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MultMesonField_hpp_
