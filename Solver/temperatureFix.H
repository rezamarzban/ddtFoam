

scalar rholow = 0.10;
int boundcounterRho = 0;

forAll(rho,i)
{
    if(rho[i]<rholow)
    {
	scalar rhofactor = max(rholow/rho[i],1.0);  // avoid negative factor
	rho[i] *= rhofactor;
	U[i]   /= rhofactor;  
	boundcounterRho++;
	
	// update k and omega
	turbulence->k()()[i] /= rhofactor;
	//turbulence->omega()[i] /= rhofactor;
	// is advected with phi, update rho*k and rho*omega not required
    }    
}

rho.correctBoundaryConditions();
U.correctBoundaryConditions();

rhoU = rho*U;

rhoTau = rho*tau;
rhofH = rho*fH;
rhoC = rho*c;

rhoE = rho*(e + 0.5*magSqr(U));

rhoEu = rho*(eu + 0.5*magSqr(U));

thermo.correct();	// should not be necessary here

p.dimensionedInternalField() = rho.dimensionedInternalField() / thermo.psi().dimensionedInternalField();
p.correctBoundaryConditions();
rho.boundaryField() = psi.boundaryField()*p.boundaryField();


reduce(boundcounterRho,sumOp<label>());  // sum values of all processors

//if(boundcounterRho) Info <<"--- increasing rho in " << boundcounterRho << " cells." << endl;

	    

