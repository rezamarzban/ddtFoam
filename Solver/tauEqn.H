
    volScalarField divrhoTauFlux = fvc::div(Riemann.rhoScalarFlux(tauIndex));


if(reconstructionStages)
{
  
    
    solve
    (
	  fvm::ddt(rho,tauHI) 
	+ divrhoTauFlux 
	- fvc::laplacian((!inviscid)*turbulence->muEff()/ScT, tau) 
	==omegaTauHI		
	,mesh.solver("tau")	// use the same solver for tau, tauHI, tauLO
    );	  
    tauHI.min(1.0);
    tauHI.max(SMALL);
    
    solve
    (
	  fvm::ddt(rho,tauLO) 
	+ divrhoTauFlux 
	- fvc::laplacian((!inviscid)*turbulence->muEff()/ScT, tau) 
	==omegaTauLO		
	,mesh.solver("tau")
    );	  
    tauLO.min(1.0);
    tauLO.max(SMALL);
    
    tau = alphaHI * tauHI + (scalar(1.0)-alphaHI) * tauLO;
    rhoTau = rho*tau;    
}
else
{
   
    solve
    (
	  fvm::ddt(rho,tau) 
	+ divrhoTauFlux 
	- fvm::laplacian((!inviscid)*turbulence->muEff()/ScT, tau) 
	==omegaTau
    );
	  
    tau.min(1.0);
    tau.max(SMALL);
    
    tauHI = tau;	// required only if the run is saved and continued with reconstructionStages>0
    tauLO = tau;
    rhoTau = rho*tau;
    
}
