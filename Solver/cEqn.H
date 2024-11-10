{

volScalarField divrhoCFlux = fvc::div(Riemann.rhoScalarFlux(cIndex));

if (!turbulentBurning) omegaC_burn *= 0.0;

omegaC = omegaC_burn + omegaC_ign;

omegaC.max(0.0);

solve
(
    fvm::ddt(rho,c) 
  + divrhoCFlux 
  - fvm::laplacian((!inviscid)*turbulence->muEff()/ScT, c) 
  ==omegaC
);

c.max(0.0);
c.min(1.0);

rhoC = rho*c;
}