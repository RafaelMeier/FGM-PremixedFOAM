/*-----------------------------------------------
Author:  Rafael Meier
Created: 11/02/2022
Update:  11/17/2022
Copyright: LabCET/UFSC

Progress Variable equation 
 
--------------------------------------------------*/

// Solve the Progress Variable equation
tmp<fvScalarMatrix> tPVEqn
    (
	(
	  fvm::ddt(rho, PV)
	+ fvm::div(phi, PV)
	  - fvm::laplacian(thermo.Dmass(), PV)  
	  == thermo.sourcePV() 
        )
    );
    
    fvScalarMatrix& PVEqn = tPVEqn.ref();

    PVEqn.relax();
    fvOptions.constrain(PVEqn);
    
    PVEqn.solve();
    fvOptions.correct(PV);
    PV = max(min( PV, 1.0 ), 0.0 );

    T = thermo.T();

    thermo.correct();
    
