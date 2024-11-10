{
// total unburned enthalpy: to get rhou and Tu

	    volScalarField divrhoHuFlux = fvc::div(Riemann.rhoEuFlux()+Riemann.pFlux());
	    
	    
	
	    solve
	    (
		fvm::ddt(rhoEu) 				
	      + divrhoHuFlux					
	      - fvc::div(sigmaDotU)				
	      - fvc::laplacian(turbulence->alphaEff(), hu) 	
	    );
	  	    
	    eu.dimensionedInternalField() = rhoEu.dimensionedInternalField() / rho.dimensionedInternalField() - 0.5*magSqr(U.dimensionedInternalField()) ;//- turbulence->k()().dimensionedInternalField();	    
            eu.correctBoundaryConditions();
            rhoEu.boundaryField() = rho.boundaryField()*
                    ( eu.boundaryField() + 0.5*magSqr(U.boundaryField()));// + turbulence->k()().boundaryField() );
		    
}