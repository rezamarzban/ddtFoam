{
// total energy equation:


          volScalarField divrhoHFlux = fvc::div(Riemann.rhoEFlux()+Riemann.pFlux());
	    
	    
	  solve
	  (
		fvm::ddt(rhoE) 			
	      + divrhoHFlux			
	      - fvc::div(sigmaDotU)		
              - fvc::laplacian(turbulence->alphaEff(), h)
	  );
	    
	    e.dimensionedInternalField() = rhoE.dimensionedInternalField() / rho.dimensionedInternalField() - 0.5*magSqr(U.dimensionedInternalField()) ;//- turbulence->k()().dimensionedInternalField();	
	    
	    if(temperatureFix)
	    {
	      	    scalar Tlow = 280.0;
		    volScalarField dummyT(T);
		    forAll(dummyT,i) dummyT[i] = Tlow;	// a volScalarField of the correct size and dimension with value = Tlow		    
		    volScalarField lowestE(e);
		    thermo.eFromT(lowestE,dummyT); // e is computed to result in a temperature of Tdummy for the local composition in each cell
		    e = max(e,lowestE);
		    rhoE = rho * (e + 0.5*magSqr(U));
	    }
	    
            e.correctBoundaryConditions();
            rhoE.boundaryField() = rho.boundaryField() * ( e.boundaryField() + 0.5*magSqr(U.boundaryField()));
		    
}