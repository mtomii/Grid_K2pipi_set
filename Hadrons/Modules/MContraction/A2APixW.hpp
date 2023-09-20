/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2APixW.hpp

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
#ifndef Hadrons_MContraction_A2APixW_hpp_
#define Hadrons_MContraction_A2APixW_hpp_

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

class A2APixWPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2APixWPar,
                                    std::string, right,
				    std::string, mes,
				    int, t_mes_base,
				    int, delt_max,
				    int, delt_min);
};

template <typename FImpl>
class TA2APixW: public Module<A2APixWPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2APixW(const std::string name);
    // destructor
    virtual ~TA2APixW(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2APixW, TA2APixW<FIMPL>, MContraction);

/******************************************************************************
 *                         TA2APixW implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2APixW<FImpl>::TA2APixW(const std::string name)
: Module<A2APixWPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2APixW<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().right, par().mes};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2APixW<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2APixW<FImpl>::setup(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  envCreate(std::vector<SpinColourVector_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2APixW<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Scalar_v;
  //typedef typename vobj::scalar_type scalar_type;
  //auto &loop  = envGet(PropagatorField, getName());
  auto &vec   = envGet(std::vector<SpinColourVector_v>, getName());
  auto &rightW = envGet(std::vector<FermionField>, par().right);
  auto &meson_f = envGet(EigenDiskVector<Complex>, par().mes);

  GridBase *grid = rightW[0].Grid();
  int nt         = env().getDim().back();
  int stepsize = par().delt_max - par().delt_min + 1;
  int v_os  = ( par().t_mes_base + par().delt_min ); //% grid->_ldimensions[orthogdim]; not necessary

  int orthogdim = 3;

  assert( grid->_ldimensions[orthogdim] >= stepsize );

  //std::vector<A2AMatrix<Scalar_s> > meson;
  std::vector<A2AMatrix<Complex> > meson;
  std::cout << "# t_mes_base: " << par().t_mes_base << std::endl;

  int tmin_rep = ( v_os  + nt ) % grid->_ldimensions[orthogdim];
  int tmax_rep = tmin_rep + stepsize - 1;
  for( int tl = tmin_rep ; tl <= tmax_rep ; ++tl ) {
    int t = tl % grid->_ldimensions[orthogdim] + grid->_lstart[orthogdim];
    int tmes;
    int t_tmes;
    for ( int t0 = 0 ; t0 < nt ; t0 += grid->_ldimensions[orthogdim] ) {
      tmes = ( t0 + par().t_mes_base ) % nt;
      t_tmes = t - tmes;
      if ( t_tmes >= par().delt_min && t_tmes <= par().delt_max ) break;
    }
    if( grid->_lstart[0] + grid->_lstart[1] + grid->_lstart[2] == 0 )
      printf("##  t_op = %d, t_mes = %d, lstart = %d\n",t,tmes,grid->_lstart[orthogdim]);
    meson.push_back(meson_f[tmes]);
  }

  int N_i = meson[0].rows();
  int N_j = meson[0].cols();
  assert ( N_j == rightW.size() );

  int rd= grid->_rdimensions[orthogdim];//2
  int e1= grid->_slice_nblock[orthogdim];//1
  int e2= grid->_slice_block [orthogdim];//64 must be 4^3
  int stride=grid->_slice_stride[orthogdim];//128
  int MFrvol = grid->_rdimensions[0]
    *          grid->_rdimensions[1]
    *          grid->_rdimensions[2]
    *          grid->_rdimensions[3]
    /          grid->_rdimensions[orthogdim]
    *          stepsize * N_i;
  vec.assign(MFrvol,Zero());
  std::cout << "Allocated " << vec.size() << " spin-color vectors for V x Pi" << std::endl;
  int volReal = e1 * e2 * stepsize * N_i;
  std::cout << "Real volume: " << e1 << " x " << e2 << " x " << stepsize << " x " << N_i << " = " << volReal << std::endl;
  std::vector<SpinColourVector_v> vec0(e1*e2*stepsize*N_i,Zero());
  //thread_for(i,N_i,{
  int i = 0;
    for(int it=0;it<stepsize;it++){
      int r = tmin_rep + it;
      int so=r*grid->_ostride[orthogdim];
      for(int j=0;j<N_j;++j){
	auto rhs_w = rightW[j].View();
	for(int n=0;n<e1;n++){
	for(int b=0;b<e2;b++){
	  int ss = so+n*stride+b;
	  //int sv = i+stepsize*(it+e1*(b+n*e2));
	  int sv = b+e2*(n+e1*(it+stepsize*i));
	  auto right = conjugate(rhs_w[ss]);
	  for(int s1=0;s1<Ns;s1++)
	  for(int c1=0;c1<Nc;c1++){
	    vec0[sv]()(s1)(c1) += right()(s1)(c1) * meson[it](i,j);
	  }
	}}
      }
    }
    //});
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2APixW_hpp_
