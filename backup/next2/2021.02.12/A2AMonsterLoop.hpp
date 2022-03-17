/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AMonsterLoop.hpp

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
#ifndef Hadrons_MContraction_A2AMonsterLoop_hpp_
#define Hadrons_MContraction_A2AMonsterLoop_hpp_

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

class A2AMonsterLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMonsterLoopPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
				    std::string, mes,
				    int, t_mes_base,
				    int, delt_max,
				    int, delt_min);
};

template <typename FImpl>
class TA2AMonsterLoop: public Module<A2AMonsterLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2AMonsterLoop(const std::string name);
    // destructor
    virtual ~TA2AMonsterLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AMonsterLoop, TA2AMonsterLoop<FIMPL>, MContraction);

/******************************************************************************
 *                     TA2AMonsterLoop implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMonsterLoop<FImpl>::TA2AMonsterLoop(const std::string name)
: Module<A2AMonsterLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AMonsterLoop<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().mes};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMonsterLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMonsterLoop<FImpl>::setup(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  //envCreateLat(PropagatorField, getName());
  envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMonsterLoop<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  //auto &loop  = envGet(PropagatorField, getName());
  auto &mat   = envGet(std::vector<SpinColourMatrix_v>, getName());
  auto &loop1 = envGet(std::vector<FermionField>, par().left);
  auto &loop2 = envGet(std::vector<FermionField>, par().right);
  auto &meson_f = envGet(EigenDiskVector<Complex>, par().mes);

  int nt         = env().getDim().back();
  int stepsize = par().delt_max - par().delt_min;
  int v_os = par().t_mes_base + par().delt_min;

  assert( nt%abs(stepsize) == 0 );

  std::vector<A2AMatrix<Complex> > meson;
  std::cout << "# t_mes_base: " << par().t_mes_base << std::endl;
  for( int t = 0 ; t < nt ; ++t ) {
    int t_mes = par().t_mes_base;
    if ( stepsize > 0 )
      t_mes = int( ((t - v_os + nt - 1)%nt)/stepsize ) * stepsize + par().t_mes_base;
    if ( stepsize < 0 )
      t_mes = ( int( (t - stepsize - v_os ) / (-stepsize) ) * (-stepsize) + par().t_mes_base ) % nt;
    std::cout << "##  t_op = " << t << ", t_mes = " << t_mes << std::endl;
    meson.push_back(meson_f[t_mes]);
  }

  GridBase *grid = loop1[0].Grid();
  int orthogdim = 3;

  int rd=grid->_rdimensions[orthogdim];//2
  int e1=    grid->_slice_nblock[orthogdim];//1
  int e2=    grid->_slice_block [orthogdim];//64 must be 4^3
  int stride=grid->_slice_stride[orthogdim];//128
  int MFrvol = grid->_rdimensions[0]
    *          grid->_rdimensions[1]
    *          grid->_rdimensions[2]
    *          grid->_rdimensions[3];
  mat.assign(MFrvol,Zero());

  int    N_i = loop1.size();
  int    N_j = loop2.size();

  std::cout << "Allocated " << mat.size() << " loop matrices" << std::endl;
  for(int i=0;i<N_i;i+=par().block)
  for(int j=0;j<N_j;j+=par().block) {
    int N_ii = MIN(N_i-i,par().block);
    int N_jj = MIN(N_j-j,par().block);
    for(int ii=0;ii<N_ii;ii+=par().cacheBlock)
    for(int jj=0;jj<N_jj;jj+=par().cacheBlock) {
      double t;
      int N_iii = MIN(N_ii-ii,par().cacheBlock);
      int N_jjj = MIN(N_jj-jj,par().cacheBlock);
      A2Autils<FImpl>::VPiW(mat, meson, &loop1[i+ii], &loop2[j+jj],
			    i+ii, j+jj, N_iii, N_jjj);
      //kernel(mat, &loop1[i+ii], &loop2[j+jj], mes, i+ii, j+jj, N_iii, N_jjj);
      /*               .
	 mat    =     / \     for each spacetime x
	              \ /
		  (a,α)x(β,b)
      */
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMonsterLoop_hpp_
