/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    fgmPremixedFoam

Description
    Flamelet generated manifold (FGM) solver for premixed combustion

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoThermo.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "lookupFGM.H"
#include "fgmThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Flamelet-generated manifold (FGM) solver for premixed combustion"
	"using density-based thermodynamics package,"
        "with enhanced buoyancy treatment."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"
  
        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        while (pimple.loop())
        {
            #include "UEqn.H"
	    #include "PVEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
               #include "pEqn.H"              
            }

        }

        rho = thermo.rho();

	Info<< "PV min/max = " << min(PV).value() << "/" << max(PV).value() << endl;
	Info<< "Source PV min/max = " << min(thermo.sourcePV()).value() << "/" << max(thermo.sourcePV()).value() << endl;
	Info<< "T min/max = " << min(T).value() << "/" << max(T).value() << endl;
	Info<< "rho min/max = " << min(rho).value() << "/" << max(rho).value() << endl;
	Info<< "p min/max = " << min(p).value() << "/ " << max(p).value() << endl;
		
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
