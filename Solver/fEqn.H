{
volScalarField divrhofHFlux = fvc::div(Riemann.rhoScalarFlux(fHIndex)); // transported scalar

	  solve
	  (
		fvm::ddt(rho,fH)
	      + divrhofHFlux 
	      - fvm::laplacian((!inviscid)*turbulence->muEff()/ScT, fH)
	  );
	  fH.max(0.0);
	  
	  rhofH = rho*fH;
}