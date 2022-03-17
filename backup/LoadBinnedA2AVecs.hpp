/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadBinnedA2AVecs.hpp

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
#ifndef Hadrons_MIO_LoadBinnedA2AVecs_hpp_
#define Hadrons_MIO_LoadBinnedA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadBinnedA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadBinnedA2AVecsPar,
                                    std::string,  filestem,
                                    bool,         multiFile,
                                    unsigned int, size,
				    unsigned int, binSize);
};

template <typename FImpl>
class TLoadBinnedA2AVecs: public Module<LoadBinnedA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadBinnedA2AVecs(const std::string name);
    // destructor
    virtual ~TLoadBinnedA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
  // read & peek
  template <typename BinnedField>
  void readpeek(std::vector<FermionField> &vec, BinnedField &bvec);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadBinnedA2AVecs, TLoadBinnedA2AVecs<FIMPL>, MIO);

/******************************************************************************
 *                     TLoadBinnedA2AVecs implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadBinnedA2AVecs<FImpl>::TLoadBinnedA2AVecs(const std::string name)
: Module<LoadBinnedA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadBinnedA2AVecs<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadBinnedA2AVecs<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadBinnedA2AVecs<FImpl>::setup(void)
{
  envCreate(std::vector<FermionField>, getName(), 1, par().size, 
  envGetGrid(FermionField));
  /*
  typedef typename FImpl::SiteSpinor::vector_type vector_type;
  envCreate(std::vector<Lattice<SiteSpinorSet<vector_type> > >, getName(), 1, par().size, 
	    envGetGrid(Lattice<SiteSpinorSet<vector_type> >));
  */
}

//template <typename FImpl, typename BinnedField>
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

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadBinnedA2AVecs<FImpl>::execute(void)
{
  auto &vec = envGet(std::vector<FermionField>, getName());

  //int bsize = NSiteSpinor;
  const int bsize = par().binSize;
  assert( vec.size() % bsize == 0 );
  int Nb = vec.size() / bsize;

  typedef typename FImpl::SiteSpinor::vector_type vector_type;

  if ( bsize == 8 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 8 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 16 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 16 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 48 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 48 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 58 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 58 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 64 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 64 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 96 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 96 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 128 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 128 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 173 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 173 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 192 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 192 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 256 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 256 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 346 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 346 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 384 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 384 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 692 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 692 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 768 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 768 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else if ( bsize == 2768 ) {
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, 2768 > SiteSpinorSet;
    std::vector<Lattice<SiteSpinorSet> >
      bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));
    readpeek(vec,bvec);
  } else {
    exit(1);
  }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadBinnedA2AVecs_hpp_
