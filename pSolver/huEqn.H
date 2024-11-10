{
// total unburned enthalpy: to get rhou and Tu
	    
	  fvScalarMatrix huEqn
	  (
	      fvm::ddt(rho,hu) 
	      + mvConvection->fvmDiv(phi, hu)
	      - fvm::laplacian(turbulence->alphaEff(), hu)
	      ==
	      DpDt*rho/rhou
	  );
	  
	  huEqn.relax();
	  huEqn.solve();  
	    
}
//uthermo.correct();	// is in the main file