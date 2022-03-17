/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AMonsterMesonField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AMonsterMesonField_hpp_
#define Hadrons_MContraction_A2AMonsterMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

BEGIN_HADRONS_NAMESPACE
//   .
//  / \
//  \ /
// --•--
//

/******************************************************************************
 *                All-to-all extended meson field creation                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AMonsterMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMonsterMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
				    int, type, // 0 for C1, C7 : 1 for C4, C10
                                    std::string, left,
                                    std::string, right,
				    std::string, loop_vw1, // new
				    std::string, loop_vw2, // new
				    std::string, mes, // external meson field
				    std::vector<int>, t_mes_os,// time of another meson field
				    int, delt_max,// signed max relative distance of 4-quark operator to meson field
                                    std::string, output,
                                    std::string, gammas1,
				    std::string, gammas2, // new
                                    std::vector<std::string>, mom);
};

class A2AMonsterMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMonsterMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma1,
				    Gamma::Algebra, gamma2);
};

template <typename T, typename FImpl>
class MonsterMesonFieldKernel: public A2AKernel<T, FImpl, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    MonsterMesonFieldKernel(const std::vector<Gamma::Algebra> &gamma1,
			    const std::vector<Gamma::Algebra> &gamma2,
			    const std::vector<LatticeComplex> &mom,
			    const int &type, const std::vector<int> &t_mes_os,
			    const int &delt_max,
			    GridBase *grid)
      // t_mes_os_(t_mes_os),
      : gamma1_(gamma1), gamma2_(gamma2), mom_(mom), type_(type), t_mes_os_(t_mes_os), delt_max_(delt_max), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
	Gt_ = grid_->GlobalDimensions()[3];
    }

    virtual ~MonsterMesonFieldKernel(void) = default;
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
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  // Quark Loop
  virtual void operator()(std::vector<Singlet_v> &cf,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      exit(1);
    }
  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      exit(1);
    }
  virtual void operator()(std::vector<ColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      exit(1);
    }
  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const int Nblock)
    {
      exit(1);
    }

  // Extended Meson Field
  virtual void operator()(A2AMatrixSet<T> &m,
			  const FieldMatrix<SpinColourVector_v > &left_loop,
			  const FermionField *right, const int offset,
			  const unsigned int orthogDim, double &t)
  {
    A2Autils<FImpl>::ExtendedMesonField(m, left_loop, right, offset, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<Singlet_v > &loop_s,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    //exit(1);
    assert( type_ == 0 );
    A2Autils<FImpl>::ExtendedMesonField(m, loop_s, left, right, gamma1_, mom_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinColourMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    //assert( type_ == 1 );
    A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, gamma2_, mom_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<ColourMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 2 );
    A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, mom_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 3 );
    A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, gamma2_, mom_, &t);
  }

  // Monstor Loop
  virtual void operator()(std::vector<Singlet_v> &cf,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 0 );
      //A2Autils<FImpl>::VPiW(cf, mes, v, w, gamma2_, i, j, Ni, Nj);
    }
  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 1 );
      //A2Autils<FImpl>::VPiW(m, mes, v, w, i, j, Ni, Nj);
    }
  virtual void operator()(std::vector<ColourMatrix_v> &cmat,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 2 );
      //A2Autils<FImpl>::VPiW(cmat, mes, v, w, gamma2_, i, j, Ni, Nj);
    }
  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 3 );
      //A2Autils<FImpl>::VPiW(m, mes, v, w, i, j, Ni, Nj);
    }

  // Left * Loop
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<SpinColourMatrix_v > &loop_m,
			  const FermionField *left, const int Nblock,
			  const int offset)
    {
      assert( type_ == 1 );
      A2Autils<FImpl>::LeftLoop(left_loop, loop_m, left, gamma1_, gamma2_, Nblock, offset);
    }
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<ColourMatrix_v > &loop_m,
			  const FermionField *left, const int Nblock,
			  const int offset)
    {
      assert( type_ == 2 );
      A2Autils<FImpl>::LeftLoop(left_loop, loop_m, left, gamma2_, Nblock, offset);
    }
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<SpinMatrix_v > &loop_m,
			  const FermionField *left, const int Nblock,
			  const int offset)
    {
      assert( type_ == 3 );
      A2Autils<FImpl>::LeftLoop(left_loop, loop_m, left, gamma1_, gamma2_, Nblock, offset);
    }

  virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
      double kernel_main;
      if ( type_ == 0 )
	kernel_main = vol_*12.0*blockSizej * 8.0 * blockSizei
	  + vol_*blockSizei*blockSizej
	  *(16.0*(6.0*3.0 + 2.0*2.0) + gamma1_.size()*(2.0*3.0+8.0));
      if ( type_ == 1 )
	kernel_main = vol_*12.0*blockSizej * 8.0 * (blockSizei+12.0)
	  + vol_*blockSizei*blockSizej * gamma1_.size()*94.0;
      return kernel_main;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
      //return Gt_*sizeof(T)*blockSizei*blockSizej*mom_.size()*gamma1_.size();
      return vol_*(12.0*sizeof(T))*blockSizei*blockSizej
	+  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej*gamma1_.size();
    }
private:
  const std::vector<Gamma::Algebra> &gamma1_;
  const std::vector<Gamma::Algebra> &gamma2_;
  const std::vector<LatticeComplex> &mom_;
  const int &type_;
  const std::vector<int> &t_mes_os_;
  const int &delt_max_;
  //const int &t_mes_os_;
  GridBase                          *grid_;
  double                            vol_;
  int Gt_;
};

template <typename FImpl>
class TA2AMonsterMesonField : public Module<A2AMonsterMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex,
				      FImpl,
                                      FermionField, 
                                      A2AMonsterMesonFieldMetadata, 
                                      HADRONS_A2AM_IO_TYPE> Computation;
  typedef MonsterMesonFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2AMonsterMesonField(const std::string name);
    // destructor
    virtual ~TA2AMonsterMesonField(void){};
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
    std::vector<Gamma::Algebra>        gamma2_;
    std::vector<std::vector<Real>>     mom_;
  int type_;
  std::vector<int> t_mes_os_;
  int delt_max_;
};

MODULE_REGISTER(A2AMonsterMesonField, ARG(TA2AMonsterMesonField<FIMPL>), MContraction);

/******************************************************************************
*               TA2AMonsterMesonField implementation                         *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMonsterMesonField<FImpl>::TA2AMonsterMesonField(const std::string name)
: Module<A2AMonsterMesonFieldPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
// Modification needed??
template <typename FImpl>
std::vector<std::string> TA2AMonsterMesonField<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().loop_vw1, par().loop_vw2, par().mes};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMonsterMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMonsterMesonField<FImpl>::setup(void)
{
    gamma1_.clear();
    gamma2_.clear();
    mom_.clear();
    if (par().gammas1 == "all")
    {
      std::cout << "gammas1 == \"all\": Not prepared yet" << std::endl;
      exit(1);
      if (par().gammas2 != "all" ) exit(1);
        gamma1_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,    
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
      if (par().gammas2 == "all") exit(1);
      gamma1_ = strToVec<Gamma::Algebra>(par().gammas1);
      gamma2_ = strToVec<Gamma::Algebra>(par().gammas2);
    }
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
    type_  = par().type;
    t_mes_os_ = par().t_mes_os;
    delt_max_ = par().delt_max;
    envCache(std::vector<ComplexField>, momphName_, 1, 
             par().mom.size(), envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
    auto &left  = envGet(std::vector<FermionField>, par().left);
    GridBase *grid = left[0].Grid();
    //envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
    envTmp(Computation, "computation", 1, grid, 
           env().getNd() - 1, mom_.size(), gamma1_.size(), par().block, 
           par().cacheBlock, this);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMonsterMesonField<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &loop1 = envGet(std::vector<FermionField>, par().loop_vw1);
    auto &loop2 = envGet(std::vector<FermionField>, par().loop_vw2);
    auto &meson_f = envGet(EigenDiskVector<Complex>, par().mes);

    //int t_mes_os      = par().t_mes_os;
    int nt         = env().getDim().back();
    assert( nt%abs(delt_max_) == 0 );
    int N_i        = left.size();
    int N_j        = right.size();
    int ngamma     = gamma1_.size();
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all MONSTER meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears1:" << std::endl;
    for (auto &g: gamma1_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Spin bilinears2:" << std::endl;
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

    //Kernel      kernel(gamma1_, ph, envGetGrid(FermionField));
    // typedef MonsterMesonFieldKernel<Complex, FImpl> Kernel;
    //Kernel      kernel(gamma1_, gamma2_, ph, type_, envGetGrid(FermionField));
    Kernel      kernel(gamma1_, gamma2_, ph, type_, t_mes_os_, delt_max_, envGetGrid(FermionField));

    auto ionameFn = [this](const unsigned int m, const unsigned int g, const unsigned int it)
    {
        std::stringstream ss;

        ss << "sort" << type_ << "_" << gamma1_[g] << "_" << gamma2_[g] << "_";
	ss << par().mes << "t" << t_mes_os_[it];
	/*
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
        {
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        }
	*/
        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const unsigned int m, const unsigned int g, const unsigned int it)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
	  + "/" + ionameFn(m, g, it) + ".h5";
    };

    //auto metadataFn = [this](const unsigned int m, const unsigned int g, const unsigned int it)
    auto metadataFn = [this](const unsigned int m, const unsigned int g)
    {
        A2AMonsterMesonFieldMetadata md;

        for (auto pmu: mom_[m])
        {
            md.momentum.push_back(pmu);
        }
        md.gamma1 = gamma1_[g];
        md.gamma2 = gamma2_[g];
	//md.t_mes_os  = t_mes_os_[it];

        return md;
    };


    for( int it = 0; it < t_mes_os_.size(); ++it ) {
      //const A2AMatrix<Complex> &meson = meson_f[t_mes_os_[it]];//const A2AMatrix<Complex>
      std::vector<A2AMatrix<Complex> > meson;
      std::cout << "# t_mes offset: " << t_mes_os_[it] << std::endl;
      for( int t = 0 ; t < nt ; ++t ) {
	int t_mes = t_mes_os_[it];
	if ( delt_max_ > 0 )
	  t_mes = ( int( (t + delt_max_ - t_mes_os_[it]) / delt_max_ ) * delt_max_ + t_mes_os_[it] ) % nt;
	if ( delt_max_ < 0 )
	  t_mes = ( int( ((t - t_mes_os_[it] + nt - 1)%nt)/(-delt_max_)) * (-delt_max_) + t_mes_os_[it] ) %nt;
	std::cout << "##  t_op = " << t << ", t_mes = " << t_mes << std::endl;
	meson.push_back(meson_f[t_mes]);
      }
      envGetTmp(Computation, computation);
      if ( type_ == 0 ) computation.execute0(left, right, loop1, loop2, meson, it, kernel, ionameFn, filenameFn, metadataFn, ngamma);
      if ( type_ == 1 ) computation.execute1(left, right, loop1, loop2, meson, it, kernel, ionameFn, filenameFn, metadataFn, ngamma);
      if ( type_ == 2 ) computation.execute2(left, right, loop1, loop2, meson, it, kernel, ionameFn, filenameFn, metadataFn, ngamma);
      if ( type_ == 3 ) computation.execute3(left, right, loop1, loop2, meson, it, kernel, ionameFn, filenameFn, metadataFn, ngamma);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMonsterMesonField_hpp_
