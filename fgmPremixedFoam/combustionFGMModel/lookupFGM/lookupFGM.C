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
       
#include "lookupFGM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lookupFGM::lookupFGM
(
 const fvMesh& mesh
)
  :
  IOdictionary
  (
    IOobject
    (
     "fgmProperties",
     mesh.time().constant(),
     mesh,
     IOobject::MUST_READ_IF_MODIFIED,
     IOobject::NO_WRITE
     )
  ),
  mesh_(mesh),
  PV_table(lookup("PV")),
  sourcePV_table(lookup("sourcePV") ),
  T_table(lookup("T") ),
  rho_table(lookup("rho"))
{

  Info << "\nFGM initialization" << endl;
  Info << "Table length: " << PV_table.size() << endl;
  Info << "Table contents:" << endl;
  Info << "{" << endl;
  Info << "PV" << endl;
  Info << "Source PV" << endl;
  Info << "Temperature" << endl;
  Info << "density" << endl;
  Info << "}\n" << endl;
 
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

lookupFGM::~lookupFGM()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::lookupFGM::interpolateValue1D
(
    const List<scalar>& table,
    scalar pvValue,
    const List<scalar>& pvTable
) const
{

    scalar interpolatedValue;

    scalar lower_pvTable=0;
    scalar upper_pvTable=0;
    scalar lower_table=0;
    scalar upper_table=0;
    
     if(pvValue == 0)
       {
	 interpolatedValue = table[0];
       }
     else
       {
	 //- A small number to prevent divide by zero
         scalar smallValue(1e-5);
         scalar rate;
 
        // INTERPOLATION ALGORITHM
        for(int j=0; j < pvTable.size(); j++ )
         {
           pvValue = min(pvValue,1.);
	  
        if(pvTable[j] >= pvValue)
	  {
	
  	    lower_pvTable = pvTable[j-1];
            upper_pvTable = pvTable[j];

	    lower_table = table[j-1];
            upper_table = table[j];
	   
	    break;
          }
	
         }
	
	 rate = (upper_table - lower_table )/ \
         max((upper_pvTable - lower_pvTable),smallValue);
  
         interpolatedValue = ( pvValue - lower_pvTable )*rate + lower_table;

         interpolatedValue = max(interpolatedValue,min(table));

        }
   
    return interpolatedValue;

}
