/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadDilutedNoise.hpp

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
#ifndef Hadrons_MIO_LoadDilutedNoise_hpp_
#define Hadrons_MIO_LoadDilutedNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadDilutedNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadDilutedNoisePar,
                                    std::string,  filestem,
                                    bool,         multiFile,
                                    unsigned int, size);
};

template <typename FImpl>
class TLoadDilutedNoise: public Module<LoadDilutedNoisePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadDilutedNoise(const std::string name);
    // destructor
    virtual ~TLoadDilutedNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadDilutedNoise, TLoadDilutedNoise<FIMPL>, MIO);

/******************************************************************************
 *                      TLoadDilutedNoise implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadDilutedNoise<FImpl>::TLoadDilutedNoise(const std::string name)
: Module<LoadDilutedNoisePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadDilutedNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadDilutedNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDilutedNoise<FImpl>::setup(void)
{
  //envCreate(std::vector<FermionField>, getName(), 1, par().size, 
  //          envGetGrid(FermionField));
  envCreate(DiluteNoise<FImpl>, getName(), 1, par().size, 
	    envGetGrid(FImpl));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDilutedNoise<FImpl>::execute(void)
{
    auto &vec = envGet(std::vector<FermionField>, getName());

    DilutedNoiseIo::read(vec, par().filestem, par().multiFile, vm().getTrajectory());
    //copyNoise(vec);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadDilutedNoise_hpp_
