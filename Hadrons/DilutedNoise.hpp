/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/DilutedNoise.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
Author: Vera Guelpers <vmg1n14@soton.ac.uk>

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
#ifndef Hadrons_DilutedNoise_hpp_
#define Hadrons_DilutedNoise_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Abstract container for diluted noise                     *
 ******************************************************************************/
template <typename FImpl>
class DilutedNoise
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    // constructor/destructor
    DilutedNoise(GridCartesian *g);
    DilutedNoise(GridCartesian *g, const unsigned int nNoise);
    virtual ~DilutedNoise(void) = default;
    // access
    std::vector<FermionField> &       getNoise(void);
    const std::vector<FermionField> & getNoise(void) const;
    const FermionField &              operator[](const unsigned int i) const;
    FermionField &                    operator[](const unsigned int i);
    void                              normalise(Real norm);
    void                              resize(const unsigned int nNoise);
    unsigned int                      size(void) const;
    GridCartesian                     *getGrid(void) const;
    // generate noise (pure virtual)
    virtual void generateNoise(GridParallelRNG &rng) = 0;
  void copyNoise(const std::vector<FermionField> &in, const int init = 0 );
private:
    std::vector<FermionField> noise_;
    GridCartesian             *grid_;
    unsigned int              nNoise_;
};

class DilutedNoiseIo
{
public:
    struct Record: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
                                        unsigned int, index);
        Record(void): index(0) {}
    };
public:
  //template <typename Field>
    //static void write(const std::string fileStem, std::vector<Field> &vec, 
    //                  const bool multiFile, const int trajectory = -1);
    template <typename Field>
    static void read(std::vector<Field> &noise, const std::string fileStem,
                     const bool multiFile, const int trajectory = -1);
private:
    static inline std::string vecFilename(const std::string stem, const int traj, 
                                          const bool multiFile)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (multiFile)
        {
            return stem + t;
        }
        else
        {
            return stem + t + ".bin";
        }
    }
};

template <typename FImpl>
class TimeDilutedSpinColorDiagonalNoise: public DilutedNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    // constructor/destructor
    TimeDilutedSpinColorDiagonalNoise(GridCartesian *g);
    virtual ~TimeDilutedSpinColorDiagonalNoise(void) = default;
    // generate noise
    virtual void generateNoise(GridParallelRNG &rng);
private:
    unsigned int nt_;
};

template <typename FImpl>
class FullVolumeSpinColorDiagonalNoise: public DilutedNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    // constructor/destructor
    FullVolumeSpinColorDiagonalNoise(GridCartesian *g, unsigned int n_src);
    virtual ~FullVolumeSpinColorDiagonalNoise(void) = default;
    // generate noise
    virtual void generateNoise(GridParallelRNG &rng);
private:
    unsigned int nSrc_;
};

template <typename FImpl>
class SparseSpinColorDiagonalNoise: public DilutedNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    // constructor/destructor
    SparseSpinColorDiagonalNoise(GridCartesian *g, unsigned int n_src, unsigned int n_sparse);
    virtual ~SparseSpinColorDiagonalNoise(void) = default;
    // generate noise
    virtual void generateNoise(GridParallelRNG &rng);
private:
    unsigned int nSrc_;
    unsigned int nSparse_;
};

/******************************************************************************
 *                    DilutedNoise template implementation                    *
 ******************************************************************************/
template <typename FImpl>
DilutedNoise<FImpl>::DilutedNoise(GridCartesian *g)
: grid_(g)
{}

template <typename FImpl>
DilutedNoise<FImpl>::DilutedNoise(GridCartesian *g,
                                  const unsigned int nNoise)
: DilutedNoise(g)
{
    resize(nNoise);
}

template <typename FImpl>
std::vector<typename DilutedNoise<FImpl>::FermionField> & DilutedNoise<FImpl>::
getNoise(void)
{
    return noise_;
}

template <typename FImpl>
const std::vector<typename DilutedNoise<FImpl>::FermionField> & DilutedNoise<FImpl>::
getNoise(void) const
{
    return noise_;
}

template <typename FImpl>
const typename DilutedNoise<FImpl>::FermionField & 
DilutedNoise<FImpl>::operator[](const unsigned int i) const
{
    return noise_[i];
}

template <typename FImpl>
typename DilutedNoise<FImpl>::FermionField & 
DilutedNoise<FImpl>::operator[](const unsigned int i)
{
    return noise_[i];
}

template <typename FImpl>
void DilutedNoise<FImpl>::normalise(Real norm)
{
    for(int i=0;i<noise_.size();i++)
    {
        noise_[i] = norm*noise_[i];
    }
}

template <typename FImpl>
void DilutedNoise<FImpl>::resize(const unsigned int nNoise)
{
    nNoise_ = nNoise;
    noise_.resize(nNoise, grid_);
}

template <typename FImpl>
unsigned int DilutedNoise<FImpl>::size(void) const
{  
    return noise_.size();
}

template <typename FImpl>
GridCartesian * DilutedNoise<FImpl>::getGrid(void) const
{
    return grid_;
}

/******************************************************************************
 *        TimeDilutedSpinColorDiagonalNoise template implementation           *
 ******************************************************************************/
template <typename FImpl>
void DilutedNoise<FImpl>::copyNoise(const std::vector<typename FImpl::FermionField> &in, const int init )
{
  std::cout << in.size() << " " << noise_.size() << std::endl;
  //assert( in.size() == noise_.size() );
  resize(in.size()-init);
  for(int i=0;i+init<in.size();++i){
    noise_[i]=in[i+init];
  }
}

/******************************************************************************
 *        TimeDilutedSpinColorDiagonalNoise template implementation           *
 ******************************************************************************/
template <typename FImpl>
TimeDilutedSpinColorDiagonalNoise<FImpl>::
TimeDilutedSpinColorDiagonalNoise(GridCartesian *g)
: DilutedNoise<FImpl>(g)
{
    nt_ = this->getGrid()->GlobalDimensions().size();
    this->resize(nt_*Ns*FImpl::Dimension);
}

template <typename FImpl>
void TimeDilutedSpinColorDiagonalNoise<FImpl>::generateNoise(GridParallelRNG &rng)
{
    typedef decltype(peekColour((*this)[0], 0)) SpinField;

    auto                       &noise = *this;
    auto                       g      = this->getGrid();
    auto                       nd     = g->GlobalDimensions().size();
    auto                       nc     = FImpl::Dimension;
    Complex                    shift(1., 1.);
    Lattice<iScalar<vInteger>> tLat(g);
    LatticeComplex             eta(g), etaCut(g);
    SpinField                  etas(g);
    unsigned int               i = 0;

    LatticeCoordinate(tLat, nd - 1);
    bernoulli(rng, eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    for (unsigned int t = 0; t < nt_; ++t)
    {
        etaCut = where((tLat == t), eta, 0.*eta);
        for (unsigned int s = 0; s < Ns; ++s)
        {
	    etas = Zero();
	    pokeSpin(etas, etaCut, s);
            for (unsigned int c = 0; c < nc; ++c)
            {
  	        noise[i] = Zero();
                pokeColour(noise[i], etas, c);
                i++;
            }
        }
    }
}

/******************************************************************************
 *        FullVolumeSpinColorDiagonalNoise template implementation           *
 ******************************************************************************/
template <typename FImpl>
FullVolumeSpinColorDiagonalNoise<FImpl>::
FullVolumeSpinColorDiagonalNoise(GridCartesian *g, unsigned int nSrc)
: DilutedNoise<FImpl>(g, nSrc*Ns*FImpl::Dimension), nSrc_(nSrc)
{}

template <typename FImpl>
void FullVolumeSpinColorDiagonalNoise<FImpl>::generateNoise(GridParallelRNG &rng)
{
    typedef decltype(peekColour((*this)[0], 0)) SpinField;

    auto                       &noise = *this;
    auto                       g      = this->getGrid();
    auto                       nd     = g->GlobalDimensions().size();
    auto                       nc     = FImpl::Dimension;
    Complex                    shift(1., 1.);
    LatticeComplex             eta(g);
    SpinField                  etas(g);
    unsigned int               i = 0;

    bernoulli(rng, eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    for (unsigned int n = 0; n < nSrc_; ++n)
    {
        for (unsigned int s = 0; s < Ns; ++s)
        {
  	    etas = Zero();
            pokeSpin(etas, eta, s);
            for (unsigned int c = 0; c < nc; ++c)
            {
	        noise[i] = Zero();
                pokeColour(noise[i], etas, c);
                i++;
            }
        }
    }
}

/******************************************************************************
 *        SparseSpinColorDiagonalNoise template implementation           *
 ******************************************************************************/
template <typename FImpl>
SparseSpinColorDiagonalNoise<FImpl>::
SparseSpinColorDiagonalNoise(GridCartesian *g, unsigned int nSrc, unsigned int nSparse)
: DilutedNoise<FImpl>(g, nSrc*Ns*FImpl::Dimension), nSrc_(nSrc), nSparse_(nSparse)
{}

template <typename FImpl>
void SparseSpinColorDiagonalNoise<FImpl>::generateNoise(GridParallelRNG &rng)
{
    typedef decltype(peekColour((*this)[0], 0)) SpinField;

    auto                       &noise = *this;
    auto                       g      = this->getGrid();
    auto                       nd     = g->GlobalDimensions().size();
    auto                       nc     = FImpl::Dimension;
    LatticeInteger             coor(g), coorTot(g); coorTot = 0.;
    Complex                    shift(1., 1.);
    LatticeComplex             eta(g), etaSparse(g);
    SpinField                  etas(g);
    unsigned int               i = 0;
    unsigned int               j = 0;
    unsigned int               nSrc_ec;
    
    if(nSrc_%nSparse_==0)
    {
         nSrc_ec = nSrc_/nSparse_;
    }
    else
    {
         nSrc_ec = (nSrc_ - nSrc_%nSparse_)/nSparse_;
    }

    for (unsigned int n = 0; n < nSrc_; ++n)
    {
        bernoulli(rng, eta);
        eta = (2.*eta - shift)*(1./::sqrt(2.));

        if(nSparse_ != 1)
        { 
        assert(g->GlobalDimensions()[1]%nSparse_ == 0);
        // # 0 # 0
        // 0 # 0 #
        // # 0 # 0
        // 0 # 0 #

        coorTot = 0;

            for(unsigned int d = 0; d < nd; ++d) 
            {
                LatticeCoordinate(coor, d);
                coorTot = coorTot + coor;
            }
            coorTot = coorTot + j;
            eta = where(mod(coorTot,nSparse_), 0.*eta, eta);
            
        }
        
        for (unsigned int s = 0; s < Ns; ++s)
        {
            etas = Zero();
            pokeSpin(etas, eta, s);
            for (unsigned int c = 0; c < nc; ++c)
            {
                noise[i] = Zero();
                pokeColour(noise[i], etas, c);
                
                i++;
                
                /**/ 
            
            }
        }
        ((n+1)%nSrc_ec == 0) ? j++: 0;
    }
    Real norm = sqrt(1./nSrc_ec);
    this->normalise(norm);
}

template <typename Field>
void DilutedNoiseIo::read(std::vector<Field> &noise, const std::string fileStem, 
			  const bool multiFile, const int trajectory)
{
    Record       record;
    ScidacReader binReader;
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < noise.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Reading noise vector " << i << std::endl;
            binReader.open(fullFilename);
            binReader.readScidacFieldRecord(noise[i], record);
	    //std::cout << vec[i] << std::endl;
            binReader.close();
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
    }
    else
    {
        binReader.open(filename);
        for (unsigned int i = 0; i < noise.size(); ++i)
        {
            LOG(Message) << "Reading noise vector " << i << std::endl;
            binReader.readScidacFieldRecord(noise[i], record);
	    //std::cout << vec[i] << std::endl;
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
        binReader.close();
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_DilutedNoise_hpp_
