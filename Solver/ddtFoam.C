/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    ddtFoam

Description
    Density based, unsteady, compressible, reactive flow solver

Authors
    Florian Ettner
    Oliver Borm  All rights reserved.
    

\*---------------------------------------------------------------------------*/

// molecular weights; required for determination of xH -> sL in XiEqn.H
#define MH2  2.01594	
#define MO2  31.9988
#define MH2O 18.0153
#define MN2  28.0134


#include "fvCFD.H"
#include "eUnburnedThermo.H"
#include "eCombustionThermo.H"
#include "turbulenceModel.H"
#include "mutUWallFunctionFvPatchScalarField.H" // yPlus
#include "ePsiChemistryModel.H"
#include "chemistrySolver.H"


#include "RiemannFluxes.H"
#include "zeroGradientFvPatchFields.H"		// for rho BC
#include "fixedRhoFvPatchScalarField.H"		// for rho BC
#include "fixedRhoEFvPatchScalarField.H"	// for rhoE BC
#include "OFstream.H"

#include "interpolationLookUpTable.H" // get species etc. from interpolation tables

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  
  #include "readChemistryProperties.H"
  #include "readInterpolationTables.H"
  #include "readAdditionalProperties.H"
  
  #include "createFields.H"
  

  #include "readTimeControls.H"
  #include "aCourantNo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    
    Info<< "\nAdaptive time-stepping is switched " ;
    if(adjustTimeStep) 
    {
        Info << "on, with CFL based on ";
        Info << " U+a (velocity + speed of sound)" << endl;
    }
    else Info << "off" << endl;

    Switch secondOrder = mesh.solutionDict().subDict("Riemann").lookupOrDefault("secondOrder",true);
    
	    Info<< " T   min...max   = " << min(T).value() << " ... " << max(T).value() << endl;
	    Info<< " p   min...max   = " << min(p).value() << " ... " << max(p).value() << endl;
	    Info<< " rho min...max   = " << min(rho).value() << " ... " << max(rho).value() << endl;    
	    Info<< " U   min...max   = " << min(mag(U)).value() << " ... " << max(mag(U)).value() << endl;
	    Info<< " tau min...max   = " << min(tau).value() << " ... " << max(tau).value() << endl;

    
    Info<< "\nStarting time loop\n" << endl;
    
    scalar flametip=-1e6;
    
    while (runTime.run())
    {
      
        #include "readTimeControls.H"
        Info << endl;
        #include "aCourantNo.H"			// based on U+a
        #include "setDeltaT.H"

        runTime++;

        Info  << "Time = " << runTime.timeName() //<< " , dt = " << runTime.deltaT().value() 
	      << ", time step = " << runTime.timeIndex() <<  endl;
	      
	volScalarField rhoRatio = rho/rhou; // for the hu eqn. I need the old rho values (before update)
	volScalarField pold=p; // required for huEqn, euEqn
	
	dimensionedScalar timestep(runTime.deltaT());  // for convection step
	
	omegaC_unquenched = uthermo.rhou() * sT * mag(fvc::grad(c)); // compute the deflagrative source term before auto-ignition, apply it afterwards
	
	omegaC_burn = omegaC_unquenched;
	if(quenching)
	{
	  omegaC_burn *= G;
	  Info << "Quenching " 
	      << 100.0* (1.0-fvc::domainIntegrate(omegaC_burn).value()/max(fvc::domainIntegrate(omegaC_unquenched).value(),SMALL))
	      << " % of the turbulent reaction " << endl;
	}
	
	if(reconstructionStages)
	{
	    #include "findMinMax_fromCellAvg.H"   
	    #include "reactMinMax.H"
	}
	else
	{
	    #include "react.H"    
	}
	// The Riemann solver uses conservative variables
	// Those variables that are solved for in the primitive formulation need to be updated in the conservative formulation
	// (actually this should be done whenever a primitive field is changed, just for safety do it here once again:)
	    rhoU   = 	rho * U;
	    rhofH  =	rho * fH;
	    rhoTau =	rho * tau;
	    rhoC   =	rho * c;
	    rhoE   =	rho * (e  + 0.5*magSqr(U)); // rhoE and rhoEu are solved for in the conservative formulation anyway
	    rhoEu   =	rho * (eu + 0.5*magSqr(U)); 
	
        // compute the convective fluxes using a Riemann solver
        Riemann.update(secondOrder);
	
	// continuity	
        volScalarField divrhoFlux  = fvc::div(Riemann.rhoFlux());
	
        solve( fvm::ddt(rho) + divrhoFlux );
	

          volTensorField tauMC("tauMC", (!inviscid)*turbulence->muEff()*dev2(fvc::grad(U)().T())); // from rhoCentralFoam	
	  
   	  surfaceScalarField a_pos = ap/(ap - am);
	  surfaceScalarField a_neg = (1.0 - a_pos);
	  surfaceScalarField sigmaDotU = 	// from rhoCentralFoam, for energy eqn.; compute here using old U
	  (	
	    (
                fvc::interpolate(turbulence->muEff())*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            ) & (a_pos*U_pos + a_neg*U_neg)*(!inviscid)
	  );
	  
	  
 	  // needed for turbulence and Xi transport eq.
	  phi = linearInterpolate(rhoU) & mesh.Sf(); // use old rhoU for simultaneous solver

	  
	// momentum
	  #include "UEqn.H"
	// mixture fraction
	  #include "fEqn.H"	  
	// auto-ignition precursor
	  #include "tauEqn.H"	  
	// reaction progress
	  #include "cEqn.H" 	  
        // turbulent flame wrinkling
	  #include "XiEqn.H"   
	// update gas composition
	  #include "lookupc.H" 
	// energy
	  #include "eEqn.H"  // updates e and h
 
	  
            // correct: thermophysical properties with new e and Y's:
	    // get new T from e, all other properties from T
	    // p does not play a role for perfect gases	and can be determined afterwards   
	    thermo.correct();
	    
            p.dimensionedInternalField() = rho.dimensionedInternalField() / thermo.psi().dimensionedInternalField();
            p.correctBoundaryConditions();
	    rho.correctBoundaryConditions();
	    
	    h = e + scalar(1.0)/psi;

	    #include "euEqn.H"	// proper BC for rhoEu required (fixed temperature on wall)
	    uthermo.correct();

	    
	    if(temperatureFix)
	    {
	       #include "temperatureFix.H"
	    }
	    	    
	    h = e + scalar(1.0)/psi;  // required for the energy eqn. in the next timestep (turbulent transport of enthalpy)
	    hu = eu + scalar(1.0)/uthermo.psiu();
	    
	    volScalarField Rgas = thermo.Cp()-thermo.Cv();
	  // speed of sound
	    cSound = sqrt((thermo.Cp()/thermo.Cv())/thermo.psi());
	    
	    //Update turbulence  
	    turbulence->correct();
	    
	    runTime.write();
	    
	    forAll(mesh.C(),i)
	    {
		if(c[i]>0.5)
		{
		    if(mesh.C()[i].component(0)>flametip)
		    {
			flametip=mesh.C()[i].component(0);
		    }
		}
	    }
	    Info << "The flame is at x= " << flametip << " m." << endl;
    }


   
    Info << "Riemann solver used: " << RiemannSolver << endl;
    
    Info << "Quenching: " << quenching << endl;
    
    Info << "\nEnd of run.\n";

    return(0);
}


// ************************************************************************* //
