// ******************************************************
// 1. compute material properties and laminar flame speed
// ******************************************************

b = scalar(1) - c;

#include "compute_sL.H"

const volScalarField& rhoUnburned = uthermo.rhou();

etaUnburned = oneKGMS * (-1.08702e-04*pow(fH,2)+1.70856e-04*fH+1.59871e-05);
etaUnburned *= pow(p/pref,0.205);
volScalarField& muUnburned = etaUnburned;
nuUnburned = etaUnburned / rhoUnburned;

alphaUnburned = oneKGMS * (-5.47712e-01*pow(fH,8)+1.44848e+00*pow(fH,7)-1.59433e+00*pow(fH,6) +9.49494e-01*pow(fH,5)-3.33735e-01*pow(fH,4)+7.10672e-02*pow(fH,3)-9.01563e-03*pow(fH,2)+5.74169e-04*fH+2.59170e-05);
alphaUnburned *= pow(p/pref,0.200);
aUnburned = alphaUnburned / rhoUnburned;


  volScalarField uPrime = sqrt(2.0/3.0*turbulence->k());
  dimensionedScalar eps0("eps0",dimensionSet(0,2,-3,0,0),1e-9);
  volScalarField epsilonTurb = max(turbulence->epsilon(),eps0);
  volScalarField lInt_Peters = 0.37*pow(uPrime,3)/epsilonTurb;	// Peters (CFX)
  //UpBySl = uPrime / sLref; // use constant sL (maximum value), otherwise UpBySl becomes too big
  //volScalarField lIntBylL = lInt_Peters/aUnburned*sLref;

  // quenching:
    nu = thermo.mu()/thermo.rho();
    gCr = pow(sL,2)/aUnburned;
    epsilonCr = 15*nu*pow(gCr,2);
    volScalarField& lInt = lInt_Peters;
    volScalarField lKolmogorov = pow(pow(nu,3)/epsilonTurb,0.25); // bei CFX wird nu eingesetzt, nicht nuUnburned
    volScalarField sigma = max(0.28*Foam::log(lInt/lKolmogorov), scalar(1e-15));
    G = 0.5 * erfc(-pow(2*sigma,-0.5)*(log(epsilonCr/epsilonTurb) + sigma/2.0));
    if(!quenching)
    {
	G = scalar(1.0);
    }
    //  Info<< "G min...max   = " << min(G).value() << " ... " << max(G).value() << endl;






// compute flame wrinkling Xi

// Calculate flame normal etc.
volVectorField n = fvc::grad(b);
volScalarField mgb = mag(n);

dimensionedScalar dMgb = 1.0e-3*
    (b*c*mgb)().weightedAverage(mesh.V())
   /((b*c)().weightedAverage(mesh.V()) + SMALL)
  + dimensionedScalar("ddMgb", mgb.dimensions(), SMALL);

mgb += dMgb;

surfaceVectorField SfHat = mesh.Sf()/mesh.magSf();
surfaceVectorField nfVec = fvc::interpolate(n);
nfVec += SfHat*(fvc::snGrad(b) - (SfHat & nfVec));
nfVec /= (mag(nfVec) + dMgb);
surfaceScalarField nf = (mesh.Sf() & nfVec);
n /= mgb;


  surfaceScalarField phiSt = fvc::interpolate(rhou*sL*Xi)*nf;
  


// Calculate coefficients for Gulder's flame speed correlation
volScalarField up = uPrimeCoef*sqrt((2.0/3.0)*turbulence->k());

volScalarField epsilonTilde = pow(uPrimeCoef, 3)*turbulence->epsilon();	// this is not the turbulent dissipation epsilon!


volScalarField tauEta = sqrt(muUnburned/(rhou*epsilonTilde));

volScalarField Reta = up/(sqrt(epsilonTilde*tauEta) + dimensionedScalar("1e-8", up.dimensions(), 1e-8));

// Calculate Xi flux
surfaceScalarField phiXi =
    phiSt
  - fvc::interpolate(fvc::laplacian(turbulence->alphaEff(), b)/mgb)*nf
  + fvc::interpolate(rho)*fvc::interpolate(sL*(1.0/Xi - Xi))*nf;


// Calculate mean and turbulent strain rates
volVectorField Ut = U + sL*Xi*n;
volScalarField Xisigmat = (n & n)*fvc::div(Ut) - (n & fvc::grad(Ut) & n);
volScalarField Xisigmas =
    ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi
  + (
        (n & n)*fvc::div(sL*n)
      - (n & fvc::grad(sL*n) & n)
    )*(Xi + scalar(1))/(2*Xi);

// Calculate Xi according to the selected flame wrinkling model
if (XiModel == "fixed")
{
    // Do nothing, Xi is fixed!
}
else if (XiModel == "algebraic")
{
    // Simple algebraic model for Xi based on Gulders correlation
    // with a linear correction function to give a plausible profile for Xi
    Xi == scalar(1) +
        (scalar(1) + (2*XiShapeCoef)*(scalar(0.5) - b))
       *XiCoef*sqrt(up/(sL + sLMin))*Reta;
}
else if (XiModel == "transport")
{
    // Calculate Xi transport coefficients based on Gulders correlation
    // and DNS data for the rate of generation
    // with a linear correction function to give a plausible profile for Xi

    volScalarField XiEqStar =
        scalar(1.001) + XiCoef*sqrt(up/(sL + sLMin))*Reta;

    volScalarField XiEq =
        scalar(1.001)
      + (scalar(1) + (2*XiShapeCoef)*(scalar(0.5) - b))
       *(XiEqStar - scalar(1.001));

    volScalarField Gstar = 0.28/tauEta;
    volScalarField XiR = Gstar*XiEqStar/(XiEqStar - scalar(1));
    volScalarField XiG = XiR*(XiEq - scalar(1.001))/XiEq;

    // Solve for the flame wrinkling
    fvScalarMatrix XiEqn
    (
        fvm::ddt(rho, Xi)
      + fvm::div(phi + phiXi, Xi, "div(phiXi,Xi)")
      - fvm::Sp(fvc::div(phiXi), Xi)
     ==
        rho*XiR
      - fvm::Sp(rho*(XiR - XiG), Xi)
      - fvm::Sp
        (
            rho*max
            (
                Xisigmat - Xisigmas,
                dimensionedScalar("0", Xisigmat.dimensions(), 0)
            ),
            Xi
        )
    );

    XiEqn.relax();
    XiEqn.solve();

    // Correct boundedness of Xi
    Xi.max(1.0);
    Info<< "max(Xi) = " << max(Xi).value() << endl;
    Info<< "max(XiEq) = " << max(XiEq).value() << endl;
    
}
else
{
    FatalError
        << args.executable() << " : Unknown Xi model " << XiModel
        << abort(FatalError);
}

sT = Xi*sL;

// compute combustion progress c = 1-b

// Create b equation
fvScalarMatrix bEqn
(
    fvm::ddt(rho, b)
  + mvConvection->fvmDiv(phi, b)
  + fvm::div(phiSt, b, "div(phiSt,b)")
  - fvm::Sp(fvc::div(phiSt), b)
  - fvm::laplacian(turbulence->alphaEff(), b)
);
bEqn.solve();
b.min(1.0);
b.max(0.0);
c = scalar(1) - b;


// not exact, for information only
omegaC = rhou * sL * Xi * mag(fvc::grad(b));

if(quenching) omegaC *= G;

