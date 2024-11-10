/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    pddtFoam


Description
    Pressure based, unsteady, compressible, reactive flow solver

Authors
    Florian Ettner
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

// molecular weights
#define MH2  2.01594	
#define MO2  31.9988
#define MH2O 18.0153
#define MN2  28.0134

#include "fvCFD.H"
#include "unburnedThermo.H"
#include "bound.H" // bounding p
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "mutUWallFunctionFvPatchScalarField.H" // yPlus
#include "hPsiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"

#include "interpolationLookUpTable.H" // get species etc. from interpolation tables

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readChemistryProperties.H"
    #include "readAdditionalProperties.H"    
    #include "readGravitationalAcceleration.H"  
    #include "readInterpolationTable.H"  
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "acousticCourantNo.H"

    #include "readPISOControls.H"


    Info << "rhoU correlation with pressure influence" << endl;
    Info << "\nThe transonic option is switched ";
        if(transonic) Info << "on" << endl;
        else Info << "off" << endl;

    Info<< "\nAdaptive time-stepping is switched " ;
    if(adjustTimeStep) 
    {
        Info << "on, with CFL based on ";
        if(!acousticCFL)
            Info << "U (flow velocity)" << endl;
        else
            Info << " U+a (flow velocity + speed of sound)" << endl;
    }
    else Info << "off" << endl;


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"
        #include "acousticCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << " , time step = " << runTime.timeIndex() << endl;


        volScalarField hOld = thermo.h();

        #include "rhoEqn.H"
	#include "convectionScheme.H" // read one convection scheme for all transported scalars
            
        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
            #include "UEqn.H"

	    #include "bEqn.H"
	    #include "fEqn.H"                    

	    #include "lookupc.H"	// same as in Riemann solver

            #include "energyEqn.H"
            thermo.correct();
	    

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
                #include "pEqn.H"
                bound(p,pmin);
            }

            //#include "createPhiC.H" // for correct computation of gradient c 
            volScalarField ctemp = 1.0*c;
            volVectorField cgradtemp = fvc::grad(ctemp);
            phiC = linearInterpolate(-cgradtemp) & mesh.Sf();

        }
        thermo.correct();

        #include "huEqn.H"
	uthermo.correct(); // calculate Tu
	
        turbulence->correct();


        cSound = sqrt((thermo.Cp()/thermo.Cv())*(thermo.Cp() - thermo.Cv())*T);


	Info<< "Combustion progress = " << 100*(scalar(1) - b)().weightedAverage(mesh.V()).value() << "%" << endl;
        Info<< " T  min...max   = " << min(T).value() << " ... " << max(T).value() << endl;
        Info<< " p  min...max   = " << min(p).value() << " ... " << max(p).value() << endl;
        Info<< " U  min...max   = " << min(mag(U)).value() << " ... " << max(mag(U)).value() << endl;
        Info<< " Ma max   = " << max(mag(U)/cSound).value() << endl;	

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
