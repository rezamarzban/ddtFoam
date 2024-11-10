
Info << "Calculating face flux field phiC\n" << endl;

volScalarField ctemp = 1.0*c;

volVectorField cgradtemp = fvc::grad(ctemp);

phiC = linearInterpolate(-cgradtemp) & mesh.Sf();

