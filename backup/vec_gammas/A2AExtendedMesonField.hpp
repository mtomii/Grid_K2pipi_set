/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AExtendedMesonField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AExtendedMesonField_hpp_
#define Hadrons_MContraction_A2AExtendedMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE
//   _
//  / \
//  \ /
// --•--
//

/******************************************************************************
 *                All-to-all extended meson field creation                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AExtendedMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
				    int, type, // 0 for C13, C23 : 1 for C16, C26
                                    std::string, left,
                                    std::string, right,
				    std::string, loop_vw1, // new
				    std::string, loop_vw2, // new
                                    std::string, output,
				    std::string, gammas1, 
                                    std::vector<std::string>, gammas2, // new
                                    std::vector<std::string>, mom);
};

class A2AExtendedMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma1,
				    Gamma::Algebra, gamma2);
};

template <typename T, typename FImpl>
class ExtendedMesonFieldKernel: public A2AKernel<T, FImpl, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
  ExtendedMesonFieldKernel(const std::vector<Gamma::Algebra> &gamma1,
			   const std::vector<std::vector<Gamma::Algebra> > &gamma2,
			   const std::vector<LatticeComplex> &mom,
			   const int &type,
			   GridBase *grid)
      : gamma1_(gamma1), gamma2_(gamma2), mom_(mom), type_(type), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~ExtendedMesonFieldKernel(void) = default;

    virtual void operator()(A2AMatrixSet<T> &m,
			    const FermionField *left, const FermionField *right,
                            const unsigned int orthogDim, double &t)
    {
      exit(1);
    }

  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iColourMatrix<vector_type> ColourMatrix_v;
  typedef iSinglet<vector_type> Singlet_v;

  // Quark Loop
  virtual void operator()(std::vector<Singlet_v> &cf,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      assert( type_ == 0 );
      //A2Autils<FImpl>::QuarkLoop(cf, v, w, gamma2_, Nblock);
    }

  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      assert( type_ == 1 );
      A2Autils<FImpl>::QuarkLoop(m, v, w, Nblock);
    }

  virtual void operator()(std::vector<ColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      assert( type_ == 2 );
      //A2Autils<FImpl>::QuarkLoop(m, v, w, gamma2_, Nblock);
    }

  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      assert( type_ == 3 );
      A2Autils<FImpl>::QuarkLoop(m, v, w, Nblock);
    }

  // Extended Meson Field
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<Singlet_v > &loop_s,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 0 );
    //A2Autils<FImpl>::ExtendedMesonField(m, loop_s, left, right, gamma1_, mom_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinColourMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 1 );
    A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, gamma2_, mom_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<ColourMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 2 );
    //A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, mom_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 3 );
    //A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, gamma2_, mom_, &t);
  }

  // Monstor Loop
  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const A2AMatrix<Complex> &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      exit(1);
    }

  virtual void operator()(std::vector<Singlet_v> &cf,
			  const FermionField *v, const FermionField *w,
			  const A2AMatrix<Complex> &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      exit(1);
    }
  virtual void operator()(std::vector<ColourMatrix_v> &cmat,
			  const FermionField *v, const FermionField *w,
			  const A2AMatrix<Complex> &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      exit(1);
    }
  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const A2AMatrix<Complex> &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      exit(1);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
      double kernel_main;
      if ( type_ == 0 )
	kernel_main = vol_*blockSizei*blockSizej
	  *(16.0*(6.0*3.0 + 2.0*2.0) + gamma1_.size()*(2.0*3.0+8.0));
      if ( type_ == 1 )
	kernel_main = vol_*blockSizei*blockSizej * gamma1_.size()*96.0;
      if ( type_ == 2 )
	kernel_main = vol_*blockSizei*blockSizej * gamma1_.size()*96.0;
      if ( type_ == 3 )
	kernel_main = vol_*blockSizei*blockSizej * gamma1_.size()*128.0;
      return kernel_main;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(12.0*sizeof(T))*blockSizei*blockSizej
	  +  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej*gamma1_.size();
    }
private:
  const std::vector<Gamma::Algebra> &gamma1_;
  const std::vector<std::vector<Gamma::Algebra> > &gamma2_;
  const std::vector<LatticeComplex> &mom_;
  const int &type_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TA2AExtendedMesonField : public Module<A2AExtendedMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex,
				      FImpl,
                                      FermionField, 
                                      A2AExtendedMesonFieldMetadata, 
                                      HADRONS_A2AM_IO_TYPE> Computation;
  typedef ExtendedMesonFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2AExtendedMesonField(const std::string name);
    // destructor
    virtual ~TA2AExtendedMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  bool                               hasPhase_{false};
  std::string                        momphName_;
  std::vector<Gamma::Algebra>        gamma1_;
  std::vector<std::vector<Gamma::Algebra> > gamma2_;
  std::vector<std::vector<Real>>     mom_;
  int type_;
};

MODULE_REGISTER(A2AExtendedMesonField, ARG(TA2AExtendedMesonField<FIMPL>), MContraction);

/******************************************************************************
*               TA2AExtendedMesonField implementation                         *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AExtendedMesonField<FImpl>::TA2AExtendedMesonField(const std::string name)
: Module<A2AExtendedMesonFieldPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
// Modification needed??
template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonField<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().loop_vw1, par().loop_vw2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::setup(void)
{
    gamma1_.clear();
    gamma2_.clear();
    mom_.clear();
    gamma1_ = strToVec<Gamma::Algebra>(par().gammas1);
    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of " 
                                + std::to_string(env().getNd() - 1));
        }
        mom_.push_back(p);
    }
    int ngamma2 = 0;
    for (auto &pstr: par().gammas2)
    {
      auto gvec = strToVec<Gamma::Algebra>(pstr);;
      gamma2_.push_back(gvec);
      ngamma2 += gvec.size();
    }
    type_ = par().type;
    envCache(std::vector<ComplexField>, momphName_, 1, 
             par().mom.size(), envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
           env().getNd() - 1, mom_.size(), ngamma2, par().block, 
           par().cacheBlock, this);
    /*
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
           env().getNd() - 1, mom_.size(), gamma1_.size(), par().block, 
           par().cacheBlock, this);
    */
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &loop1 = envGet(std::vector<FermionField>, par().loop_vw1);
    auto &loop2 = envGet(std::vector<FermionField>, par().loop_vw2);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int ngamma     = gamma1_.size();
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all EXTENDED meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Gamma matrices 1:" << std::endl;
    for (auto &g: gamma1_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Gamma matrices 2:" << std::endl;
    for (auto &g: gamma2_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (unsigned int j = 0; j < nmom; ++j)
        {
            Complex           i(0.0,1.0);
            std::vector<Real> p;

            envGetTmp(ComplexField, coor);
            ph[j] = Zero();
            for(unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        hasPhase_ = true;
        stopTimer("Momentum phases");
    }

    auto ionameFn = [this](const unsigned int m, const unsigned int g0)
    {
        std::stringstream ss;
	int g2 = g0;
	int g  = 0;
	while ( g2 >= gamma2_[g].size() ) {
	  g2 -= gamma2_[g].size();
	  g++;
	}
	assert ( g < gamma1_.size() );

        ss << "type" << type_ << "_" << gamma1_[g] << "_" << gamma2_[g][g2] << "_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
        {
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        }

        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const unsigned int m, const unsigned int g0)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
	  + "/" + ionameFn(m, g0) + ".h5";
    };

    auto metadataFn = [this](const unsigned int m, const unsigned int g0)
    {
        A2AExtendedMesonFieldMetadata md;

        for (auto pmu: mom_[m])
        {
            md.momentum.push_back(pmu);
        }

	int g2 = g0;
	int g  = 0;
	while ( g2 >= gamma2_[g].size() ) {
	  g2 -= gamma2_[g].size();
	  g++;
	}
	assert ( g < gamma1_.size() );
        md.gamma1 = gamma1_[g];
        //md.gamma2 = gamma1_[g];
        md.gamma2 = gamma2_[g][g2];
        
        return md;
    };

    //Kernel      kernel(gamma1_, ph, envGetGrid(FermionField));
    // typedef ExtendedMesonFieldKernel<Complex, FImpl> Kernel;
    Kernel      kernel(gamma1_, gamma2_, ph, type_, envGetGrid(FermionField));

    envGetTmp(Computation, computation);
    if ( type_ == 0 ) computation.execute0(left, right, loop1, loop2, kernel, ionameFn, filenameFn, metadataFn, ngamma);
    if ( type_ == 1 ) computation.execute1(left, right, loop1, loop2, kernel, ionameFn, filenameFn, metadataFn);
    if ( type_ == 2 ) computation.execute2(left, right, loop1, loop2, kernel, ionameFn, filenameFn, metadataFn, ngamma);
    if ( type_ == 3 ) computation.execute3(left, right, loop1, loop2, kernel, ionameFn, filenameFn, metadataFn);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AExtendedMesonField_hpp_
