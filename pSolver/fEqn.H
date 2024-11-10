fvScalarMatrix fEqn
(
    fvm::ddt(rho, fH)
  + mvConvection->fvmDiv(phi, fH)
  - fvm::laplacian(turbulence->muEff()/ScT, fH)
);

fEqn.relax();
fEqn.solve();

