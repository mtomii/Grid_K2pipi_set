/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSolver/A2AVectorsV.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: fionnoh <fionnoh@gmail.com>

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
#ifndef Hadrons_MSolver_A2AVectorsV_hpp_
#define Hadrons_MSolver_A2AVectorsV_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsVPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsVPar,
                                  std::string, win,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver,
                                  std::string, output,
                                  bool,        multiFile);
		      //int, nevecs,
};

template <typename FImpl, typename Pack>
class TA2AVectorsV : public Module<A2AVectorsVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2AVectorsV(const std::string name);
    // destructor
    virtual ~TA2AVectorsV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
  unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(A2AVectorsV, 
    ARG(TA2AVectorsV<FIMPL, BaseFermionEigenPack<FIMPL>>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsV, 
    ARG(TA2AVectorsV<ZFIMPL, BaseFermionEigenPack<ZFIMPL>>), MSolver);

/******************************************************************************
 *                       TA2AVectorsV implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectorsV<FImpl, Pack>::TA2AVectorsV(const std::string name)
: Module<A2AVectorsVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectorsV<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;

    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver + sub_string);
    in.push_back(par().win);

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectorsV<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectorsV<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &win        = envGet(std::vector<FermionField>, par().win);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver + sub_string);
    int         Ls          = env().getObjectLs(par().action);
    //Nl_ = par().nevecs;

    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    /*
    if (hasLowModes)
      envCreateDerived(BasePack, Pack, "evecs", Ls, Nl_,
		       envGetRbGrid(Field, Ls), gridIo);
    */
    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nl_ + win.size(), envGetGrid(FermionField));
    envCreateDerived(DilutedNoise<FImpl>, 
                     TimeDilutedSpinColorDiagonalNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectorsV<FImpl, Pack>::execute(void)
{
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver + sub_string);
    auto        &win       = envGet(std::vector<FermionField>, par().win);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    //auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);
    auto        &noise     = envGet(DilutedNoise<FImpl>, getName() );
    noise.copyNoise(win);
    envGetTmp(A2A, a2a);

    //assert( Nl_ == 0 );
    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using eigenpack ("
                     << Nl_ << " low modes) and noise ("
		     << noise.size() << " noise vectors) from W vectors '"
		     << par().win << "'" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using noise '" << par().win << "' (" << noise.size() 
                     << " noise vectors)" << std::endl;
    }
    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
      auto &epack  = envGet(Pack, par().eigenPack );

        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
	  a2a.makeLowModeV(v[il], epack.evec[il], epack.eval[il]);
	  //a2a.makeLowModeV(v[il], win[il], win[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[il], f5, epack.evec[il], epack.eval[il]);
            //a2a.makeLowModeV5D(v[il], f5, epack.evec[il], win[il]);
        }
        stopTimer("V low mode");
    }

    // High modes
    for (unsigned int ih = 0; ih < noise.size(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[Nl_ + ih], noise[ih]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise[ih]);
        }
        stopTimer("V high mode");
    }

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsV_hpp_
