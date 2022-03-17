/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadQlatA2AVecs.hpp

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
#ifndef Hadrons_MIO_LoadQlatA2AVecs_hpp_
#define Hadrons_MIO_LoadQlatA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

#include <qlat/qlat.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadQlatA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadQlatA2AVecsPar,
                                    std::string,  filestem,
                                    unsigned int, size);
};

template <typename FImpl>
class TLoadQlatA2AVecs: public Module<LoadQlatA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadQlatA2AVecs(const std::string name);
    // destructor
    virtual ~TLoadQlatA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
  // read & peek
  //template <typename BinnedField>
  //void readpeek(std::vector<FermionField> &vec, BinnedField &bvec);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadQlatA2AVecs, TLoadQlatA2AVecs<FIMPL>, MIO);

/******************************************************************************
 *                      TLoadQlatA2AVecs implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadQlatA2AVecs<FImpl>::TLoadQlatA2AVecs(const std::string name)
: Module<LoadQlatA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadQlatA2AVecs<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadQlatA2AVecs<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadQlatA2AVecs<FImpl>::setup(void)
{
  envCreate(std::vector<FermionField>, getName(), 1, par().size, 
  envGetGrid(FermionField));
  /*
  typedef typename FImpl::SiteSpinor::vector_type vector_type;
  envCreate(std::vector<Lattice<SiteSpinorSet<vector_type> > >, getName(), 1, par().size, 
	    envGetGrid(Lattice<SiteSpinorSet<vector_type> >));
  */
}

#if 0
//template <typename FImpl, typename QlatField>
template<typename FImpl>
template <typename BinnedField>
void TLoadBinnedA2AVecs<FImpl>::readpeek(std::vector<FermionField> &vec, BinnedField &bvec)
{
  const int bsize = par().binSize;
  A2AVectorsIo::read(bvec, par().filestem, par().multiFile, vm().getTrajectory());
  std::cout << "Vecsize: " << vec.size() << ", BVecsize: " << bvec.size() << std::endl;
  for( int i = 0 ; i < vec.size() ; i += bsize ) {
    for( int j = 0 ; j < bsize ; ++j ) {
      vec[i+j] = peekLorentz(bvec[i/bsize],j);
    }
  }
}
#endif

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadQlatA2AVecs<FImpl>::execute(void)
{
  auto &vec = envGet(std::vector<FermionField>, getName());

  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourVector<scalar_type> SpinColourVector_s;


  GridBase *grid = vec[0].Grid();
  qlat::begin(UniqueID(),qlat::Coordinate(GJP.Xnodes(),GJP.Ynodes(),GJP.Znodes(),GJP.Tnodes()));
  qlat::Coordinate total_site(grid->_gdimensions[0],grid->_gdimensions[1],grid->_gdimensions[2],grid->_gdimensions[3]);
  qlat::Geometry geo;
  geo.init(total_site, 1);
  qlat::FieldM<qlat::WilsonVector, 1> fvec;
  fvec.init(geo);

  typedef typename FImpl::SiteSpinor::vector_type vector_type;
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadQlatA2AVecs_hpp_
