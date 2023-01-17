/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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
 
\*---------------------------------------------------------------------------*/

#include "fgmThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fgmThermo::fgmThermo
(
    const fvMesh& mesh,
    const lookupFGM& lookupFGM,
    volScalarField& p,
    volScalarField& PV
)
:

  IOdictionary
  (
    IOobject
    (
     "thermophysicalProperties",
     mesh.time().constant(),
     mesh,
     IOobject::MUST_READ_IF_MODIFIED,
     IOobject::NO_WRITE
     )
  ),

  fgmThermoModel_(lookup("thermoType")),
  Le_(lookupOrDefault<scalar>("Le",1.0)),
  Rgas_(lookupOrDefault<scalar>("Rgas",287.0)),
  
  fgmTable_(lookupFGM),

  p_(p),
  
  PV_(PV),

  T_
   (
       IOobject
       (
           "T",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionSet(0,0,0,1,0,0,0)
    ),

  rho_
   (
       IOobject
       (
           "thermo:rho",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimDensity
    ),
  
  psi_
   (
    IOobject
    (
        "thermo:psi",
        mesh.time().timeName(),
        mesh
    ),
    mesh,
    dimensionSet(0,-2,2,0,0,0,0)
    ),
  
   Dmass_
   (
    IOobject
    (
        "thermo:Dmass",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime,1E-5)
    ),

  mu_
   (
    IOobject
    (
        "thermo:mu",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime,1E-5)
    ),
  
  sourcePV_
  (
    IOobject
    (
        "sourcePV",
        mesh.time().timeName(),
        mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume/dimTime)
   ),

  dpdt_(lookupOrDefault<Switch>("dpdt", false))
  
{

  
  // Checking the FGM Model
  if (fgmThermoModel_ == "fgm1DModelDNS")
    {
       Info << "\nInitializing fgmThermo\n" << endl;
    }
  else
    {
      FatalErrorIn(
		   "fgmThermo::fgmThermoModel"
		   )
	<< "The model " << fgmThermoModel_ << " does not exist." << abort(FatalError); 

    }
  
   correct();
  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fgmThermo::~fgmThermo()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fgmThermo::correct()
{
  scalarField& TCells = T_.primitiveFieldRef();
  scalarField& psiCells = psi_.primitiveFieldRef();
  scalarField& sourcePVCells = sourcePV_.primitiveFieldRef();
  scalarField& rhoCells = rho_.primitiveFieldRef();
  scalarField& DmassCells = Dmass_.primitiveFieldRef();
  scalarField& muCells = mu_.primitiveFieldRef();
 
  const scalarField& pCells = p_.internalField();
  const scalarField& PVCells = PV_.internalField();

  // Interpolate for internal field
  forAll(TCells, celli)
  {
      TCells[celli] = fgmTable_.interpolateValue1D
                        (
			    fgmTable_.T_table,
			    PVCells[celli],
			    fgmTable_.PV_table
                        );

      DmassCells[celli] = massDiffusivity_model(TCells[celli]);

      muCells[celli] = viscosity_model(TCells[celli]);
	
      rhoCells[celli] = fgmTable_.interpolateValue1D
                        (
			    fgmTable_.rho_table,
			    PVCells[celli],
			    fgmTable_.PV_table
                        );

      psiCells[celli] = compressibility_model(rhoCells[celli], pCells[celli]);
      
      sourcePVCells[celli] = fgmTable_.interpolateValue1D
                        (
			    fgmTable_.sourcePV_table,
			    PVCells[celli],
			    fgmTable_.PV_table
                        );
    
  }

  // Interpolate for patches
  forAll(T_.boundaryField(), patchi)
    {
      fvPatchScalarField& pT = T_.boundaryFieldRef()[patchi];
      fvPatchScalarField& ppsi = psi_.boundaryFieldRef()[patchi];
      fvPatchScalarField& prho = rho_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pDmass = Dmass_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pmu = mu_.boundaryFieldRef()[patchi];
      fvPatchScalarField& psourcePV = sourcePV_.boundaryFieldRef()[patchi];

      const fvPatchScalarField& pp = p_.boundaryField()[patchi];
      const fvPatchScalarField& pPV = PV_.boundaryField()[patchi];

      forAll(pT, facei)
      {
                pT[facei] = fgmTable_.interpolateValue1D
                                (
		                    fgmTable_.T_table,
			            pPV[facei],
			            fgmTable_.PV_table
                                );

		pDmass[facei] = massDiffusivity_model(pT[facei]);

		pmu[facei] = viscosity_model(pT[facei]);

		prho[facei] = fgmTable_.interpolateValue1D
                                (
		                    fgmTable_.rho_table,
			            pPV[facei],
			            fgmTable_.PV_table
                                );

		ppsi[facei] = compressibility_model(prho[facei],pp[facei]);
		
		psourcePV[facei] = fgmTable_.interpolateValue1D
                                (
		                    fgmTable_.sourcePV_table,
			            pPV[facei],
			            fgmTable_.PV_table
                                );	
      }
    }
}

void Foam::fgmThermo::correctRho
(
 const Foam::volScalarField& deltaRho
)
{
  /*
   Density variation due to the dynamic pressure is neglected.
   Low Mach approximation
  */

  /*
  const dimensionedScalar rhoMax("rhoMax", dimDensity, max(fgmTable_.rho_table));
  const dimensionedScalar rhoMin("rhoMin", dimDensity, min(fgmTable_.rho_table));
  rho_ += deltaRho;
  rho_ = max(rho_, rhoMin);
  rho_ = min(rho_, rhoMax);
  */  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fgmThermo::compressibility_model(const scalar rho, const scalar p) const
{
  return rho/p;   
}

Foam::scalar Foam::fgmThermo::massDiffusivity_model(const scalar T) const
{

  scalar T298 = 298;
  scalar C069 = 0.69;
  scalar CD = 2.58E-5;

  return CD*pow(T/T298,C069)/Le_; 
  
}

Foam::scalar Foam::fgmThermo::viscosity_model(const scalar T) const
{
  // Model: Sutherland's law 
  scalar muRef = 1.7894E-5;
  scalar TRef = 273.15;
  scalar S = 110.4;
 
  return muRef*pow(T/TRef,1.5)*((TRef+S)/(T+S)); 
  
}

// ************************************************************************* //
