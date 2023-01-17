/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
 
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
 
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
Class
    Foam::fgmThermo
 
Description
    Class to compute all transport and thermodynamic data according
    FGM model.

    @author: Rafael Becker Meier
    @version: 11/28/2022
 
SourceFiles
    fgmThermo.C
 
\*---------------------------------------------------------------------------*/

#ifndef fgmThermo_H
#define fgmThermo_H

#include "IOdictionary.H"

#include "lookupFGM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class fgmThermo
  :
    public IOdictionary
{

  private:
  
     // Protected data

       //- FGM Model
       word fgmThermoModel_;

      //- Lewis number [-]
       const scalar Le_;
  
       //- Constant of gas mixture [J/Kg K]
       const scalar Rgas_;
  
       //- Look up table and routines
       const lookupFGM& fgmTable_;

  // Fields

       //- Pressure [Pa]
       volScalarField p_;
  
       //- Progress variable [-]
       volScalarField& PV_;

       //- Temperature [K]
       volScalarField T_;

       //- Density [Kg/m3]
       volScalarField rho_;

       //- Compressibility at constant temperature [s2/m2]
       volScalarField psi_;

       //- Mass diffusivity [Kg/m s]
       volScalarField Dmass_;

       //- Viscosity [Kg/m s]
       volScalarField mu_;
  
       //- Rate reaction of progress variable [Kg/m3s]
       volScalarField sourcePV_;

       //- Should the dpdt term be included in the enthalpy equation
       Switch dpdt_;
       
public:

        //- Return the instantiated type name
        static word typeName()
          {
            return "fgmThermo";
           }
  
    // Constructors
    fgmThermo
     (
       const fvMesh&,
       const lookupFGM&,
       volScalarField&,
       volScalarField&
     );

    //- Destructor
    virtual ~fgmThermo();
  
    // Member functions

        //- Calculate the thermo variables
        void correct();

        //- Add the given density correction to the density field.
        //  Used to update the density field following pressure solution
        void correctRho(const volScalarField&);

        //- Calculate compressibility rho/p
        scalar compressibility_model(const scalar, const scalar) const;
  
        //- Calculate the mass diffusion coefficient
        scalar massDiffusivity_model(const scalar) const;

        //- Calculate the viscosity coefficient: Sutherland's law
        scalar viscosity_model(const scalar) const;

        // Return private fields
  
        //- In this model dpdt term is always neglected
        Switch dpdt() const
        {
          return dpdt_;
        }

        //- Temperature [K]
        virtual tmp<volScalarField> T() const
	{
	  return T_;
	}
    
         // Return non-const access to the local temperature field [K]
	virtual volScalarField& T()
	{
	  return T_;
	}
  
        //- Pressure [Pa]
        virtual tmp<volScalarField> p() const
	{
	  return p_;
	}

        // Return non-const access to the local pressure field [Pa]
        virtual volScalarField& p()
	{
	  return p_;
	}

        //- Density [kg/m^3]
	virtual tmp<volScalarField> rho() const
	{
	  return rho_;
	}

        //- Return non-const access to the local density field [kg/m^3]
	virtual volScalarField& rho()
	{
	  return rho_;
	}

         //- Compressibility [s2/m2]
        virtual tmp<volScalarField> psi() const
	{
	  return psi_;
	}

        //- Return non-const access to the local compressibility field [s2/m2]
	virtual volScalarField& psi()
	{
	  return psi_;
	}
  
        //- Rate reaction of progress variable [Kg/m3s]
        virtual tmp<volScalarField> sourcePV() const
	{
	  return sourcePV_;
	}

        //- Return non-const access to the local rate reaction of
        // progress variable field [kg/m3s]
	virtual volScalarField& sourcePV()
	{
	  return sourcePV_;
	}

        //- Mass diffusivity [Kg/ms]
        virtual tmp<volScalarField> Dmass() const
	{
	  return Dmass_;
	}

        //- Return non-const access to the local mass diffusivity field [kg/ms]
        virtual volScalarField& Dmass()
	{
	  return Dmass_;
	}

        //- Dynamic viscosity [kg/m/s]
        virtual tmp<volScalarField> mu() const
	{
	  return mu_;
	}

        //- Return non-const access to the local dynamic viscosity field [kg/ms]
        virtual volScalarField& mu()
	{
	  return mu_;
	}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
