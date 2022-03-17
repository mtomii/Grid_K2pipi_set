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
				    std::string, propField,
#if 0
				    std::string, loop_vw1, // new
				    std::string, loop_vw2, // new
				    std::string, mes, // external meson field
				    std::vector<int>, t_mes_os,
				    int, delt_max,
				    int, delt_min,
#endif
                                    std::string, output,
                                    std::string, gammas1,
				    std::string, gammas2);
};

class A2AMonsterMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMonsterMesonFieldMetadata,
				    std::string, propF,
                                    std::string, gamma1,
				    std::string, gamma2);
};

template <typename T, typename FImpl>
class MonsterMesonFieldKernel: public A2AKernel<T, FImpl, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
  MonsterMesonFieldKernel(const std::vector<std::vector<Gamma::Algebra> > &gamma1,
			  const std::vector<std::vector<Gamma::Algebra> > &gamma2,
			  const std::vector<std::string> &nameg1,
			  const std::vector<std::string> &nameg2,
			  const int &type,
			  const std::string &propF,
			  GridBase *grid)
      // t_mes_os_(t_mes_os),
    : gamma1_(gamma1), gamma2_(gamma2), nameg1_(nameg1), nameg2_(nameg2), type_(type), propField_(propF), grid_(grid)
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
    A2Autils<FImpl>::ExtendedMesonField(m, loop_s, left, right, gamma1_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinColourMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 1 );
    //A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, gamma2_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<ColourMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 2 );
    //A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, &t);
  }
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinMatrix_v > &loop_m,
			  const FermionField *left, const FermionField *right,
			  const unsigned int orthogDim, double &t)
  {
    assert( type_ == 3 );
    //A2Autils<FImpl>::ExtendedMesonField(m, loop_m, left, right, gamma1_, gamma2_, &t);
  }

  // Monster Loop
  virtual void operator()(std::vector<Singlet_v> &cf,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 0 );
      A2Autils<FImpl>::VPiW(cf, mes, v, w, gamma2_, i, j, Ni, Nj);
    }
  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 1 );
      A2Autils<FImpl>::VPiW(m, mes, v, w, i, j, Ni, Nj);
    }
  virtual void operator()(std::vector<ColourMatrix_v> &cmat,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 2 );
      A2Autils<FImpl>::VPiW(cmat, mes, v, w, gamma2_, i, j, Ni, Nj);
    }
  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const FermionField *v, const FermionField *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj)
    {
      assert( type_ == 3 );
      A2Autils<FImpl>::VPiW(m, mes, v, w, i, j, Ni, Nj);
    }

  // Left * Loop
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<SpinColourMatrix_v > &loop_m,
			  const FermionField *left, const int Nblock,
			  const int offset)
    {
      if ( type_ == 0 ) A2Autils<FImpl>::LeftLoop1(left_loop, loop_m, left, gamma1_, gamma2_, Nblock, offset);
      if ( type_ == 1 ) A2Autils<FImpl>::LeftLoop1(left_loop, loop_m, left, gamma1_, gamma2_, Nblock, offset);
      if ( type_ == 2 ) A2Autils<FImpl>::LeftLoop2(left_loop, loop_m, left, gamma1_, gamma2_, Nblock, offset);
      if ( type_ == 3 ) A2Autils<FImpl>::LeftLoop3(left_loop, loop_m, left, gamma1_, gamma2_, Nblock, offset);
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
      //return Gt_*sizeof(T)*blockSizei*blockSizej*gamma1_.size();
      return vol_*(12.0*sizeof(T))*blockSizei*blockSizej
	+  vol_*(2.0*sizeof(T)*blockSizei*blockSizej*gamma1_.size());
    }
private:
  const std::vector<std::vector<Gamma::Algebra> > &gamma1_;
  const std::vector<std::vector<Gamma::Algebra> > &gamma2_;
  const std::vector<std::string> &nameg1_;
  const std::vector<std::string> &nameg2_;
  const int &type_;
  const std::string &propField_;
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
  std::vector<std::vector<Gamma::Algebra> >       gamma1_;
  std::vector<std::vector<Gamma::Algebra> >       gamma2_;
  std::vector<std::string> nameg1_;
  std::vector<std::string> nameg2_;
  int type_;
  std::string propField_;
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
  std::vector<std::string> in = {par().left, par().right, par().propField};

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
    type_  = par().type;
    propField_ = par().propField;
    envCache(std::vector<ComplexField>, momphName_, 1, 
             1, envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
    auto &left  = envGet(std::vector<FermionField>, par().left);
    GridBase *grid = left[0].Grid();
    //envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
    envTmp(Computation, "computation", 1, grid, 
           env().getNd() - 1, 1, gamma1_.size(), par().block, 
           par().cacheBlock, this);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMonsterMesonField<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;

  auto &left  = envGet(std::vector<FermionField>, par().left);
  auto &right = envGet(std::vector<FermionField>, par().right);
  auto &loop  = envGet(std::vector<SpinColourMatrix_v>, par().propField);

  int nt         = env().getDim().back();

  int N_i        = left.size();
  int N_j        = right.size();
  int block      = par().block;
  int cacheBlock = par().cacheBlock;

  int ngamma = 0;
  for ( int ig = 0 ; ig < gamma1_.size() ; ++ig ) {
    ngamma += gamma1_[ig].size();
  }

  LOG(Message) << "Computing all-to-all MONSTER meson fields" << std::endl;
  LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
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
	       << "/bilinear)" << std::endl;

  //auto &ph = envGet(std::vector<ComplexField>, momphName_);

  Kernel kernel(gamma1_, gamma2_, nameg1_, nameg2_, type_, propField_, envGetGrid(FermionField));

  auto ionameFn = [this](const unsigned int g)
  {
    std::stringstream ss;

    ss << "sort" << type_ << "_" << nameg1_[g] << "_" << nameg2_[g] << "_";
    ss << propField_;
    return ss.str();
  };

  auto filenameFn = [this, &ionameFn](const unsigned int g)
  {
    return par().output + "." + std::to_string(vm().getTrajectory()) 
      + "/" + ionameFn(g) + ".h5";
  };

  auto metadataFn = [this](const unsigned int g)
  {
    A2AMonsterMesonFieldMetadata md;
    md.propF = propField_;
    md.gamma1 = nameg1_[g];
    md.gamma2 = nameg2_[g];

    return md;
  };

  envGetTmp(Computation, computation);
  computation.execute1(left, right, loop, kernel, ionameFn, filenameFn, metadataFn);
#if 0
  if ( type_ == 0 ) computation.execute0(left, right, loop, kernel, ionameFn, filenameFn, metadataFn, ngamma);
  if ( type_ == 1 ) computation.execute1(left, right, loop, kernel, ionameFn, filenameFn, metadataFn);
  if ( type_ == 2 ) computation.execute2(left, right, loop1, loop2, meson, it, kernel, ionameFn, filenameFn, metadataFn, ngamma);
  if ( type_ == 3 ) computation.execute3(left, right, loop1, loop2, meson, it, kernel, ionameFn, filenameFn, metadataFn);
#endif
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMonsterMesonField_hpp_
