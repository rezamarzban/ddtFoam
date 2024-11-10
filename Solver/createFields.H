
Info<< "createFields:" << endl;

const dimensionedScalar oneMS("oneMS", dimensionSet(0,1,-1,0,0), 1.0); // for sL correlation
const dimensionedScalar oneM("oneM", dimensionSet(0,1,0,0,0), 1.0); // for Zbikowski correlation
const dimensionedScalar oneKGM3("oneKGM3", dimensionSet(1,-3,0,0,0), 1.0); // for rhoUnburned
const dimensionedScalar oneM2S("oneM2S", dimensionSet(0,2,-1,0,0), 1.0); // for aUnburned / nuUnburned
const dimensionedScalar oneKGMS("oneKGMS", dimensionSet(1,-1,-1,0,0), 1.0); // for alphaUnburned / muUnburned
const dimensionedScalar pref("pref", dimensionSet(1,-1,-2,0,0), 1.0e5); // for rhoUnburned/alphaUnburned/muUnburned
const dimensionedScalar pmin("pmin", dimensionSet(1,-1,-2,0,0), 0.2e5); // for bounding
const dimensionedScalar T0("T0", dimensionSet(0,0,0,1,0), 293.0); // for calculation of Tunburned
const dimensionedScalar p0("p0", dimensionSet(1,-1,-2,0,0), 1.0e5); // for calculation of Tunburned

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<ePsiChemistryModel> pChemistry
(
    ePsiChemistryModel::New(mesh)
);

ePsiChemistryModel& chemistry = pChemistry();
eCombustionThermo& thermo = chemistry.thermo();

basicMultiComponentMixture& composition = thermo.composition();

Info << "thermo: " << thermo << endl; 
Info << "species: " << composition.species() << endl; 

PtrList<volScalarField>& Y = composition.Y();

    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

word inertSpecie(thermo.lookup("inertSpecie"));

#include "setupInertSpecie.H"
    volScalarField& e = thermo.e();		
    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();			
    const volScalarField& T = thermo.T();

    
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    e+scalar(1.0)/psi
);
//h.write();




    

#include "rhoBoundaryTypes.H"
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho(),
	rhoBoundaryTypes
    );
    rho.write();
    

    //Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    U.write();

#   include "compressibleCreatePhi.H" 
    //phi.write();

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );


    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );
    //rhoU.write();
    
volScalarField fH
(
    IOobject
    (
        "fH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField rhofH
(
    IOobject
    (
        "rhofH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*fH	// results in a calculated BC
);
//rhofH.write();
    
#include "rhoEBoundaryTypes.H"
    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(e + 0.5*magSqr(U))//,// + turbulence->k()),
	//rhoEBoundaryTypes	// wäre nicht unbedingt nötig, da RB-Werte aus rho,e,U-RB berechnet werden können
    );
    //rhoE.write();
    //e.write();
    
    
    eUnburnedThermo uthermo(mesh, p, fH, yO2inAir);
    const volScalarField& Tu=uthermo.Tu();
    //Info << " Tu varies between " << min(Tu).value() << " and " << max(Tu).value() << " K." << endl;
    const volScalarField& rhou = uthermo.rhou();
    volScalarField& eu = uthermo.eu();		
    volScalarField hu = eu + scalar(1.0) / uthermo.psiu();
    volScalarField rhoEu
    (
        IOobject
        (
            "rhoEu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho*(eu + 0.5*magSqr(U))//,// + turbulence->k()),
    );
    rhoEu.write();

    
    
    
    


label h2Index=-1;
label o2Index=-1;
label n2Index=-1;
forAll(Y, i)
{
	if(Y[i].name()=="H2") h2Index=i;	
	if(Y[i].name()=="O2") o2Index=i;
	if(Y[i].name()=="N2") n2Index=i;
}
// alternative:
// h2Index = composition.species()["H2"];



volScalarField tau  // auto-ignition process
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
tau.max(0.0);
tau.min(1.0);
tau.write();

volScalarField rhoTau
(
    IOobject
    (
        "rhoTau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*tau
    //,rhoBoundaryTypes // vorläufig
);


volScalarField omegaTau  // auto-ignition source term
(
    IOobject
    (
        "omegaTau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaTau", dimMass/dimVolume/dimTime, 0.0)
);





volScalarField tIgn  // auto-ignition delay
(
    IOobject
    (
        "tIgn",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tIgn_startValue", dimTime, 1e12),
    zeroGradientFvPatchScalarField::typeName // use zeroGradient instead of calculated BC
);


volScalarField c  // TFC reaction progress
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

c.max(0.0);
c.min(1.0);



volScalarField rhoC
(
    IOobject
    (
        "rhoC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*c
);



    volScalarField Xi
    (
        IOobject
        (
            "Xi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



dimensionedScalar sLMin("sLMin",dimensionSet(0,1,-1,0,0),0.001);    

volScalarField sL  // laminar flame speed
(
    IOobject
    (
        "sL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sL", dimLength/dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName // use zeroGradient instead of calculated BC 
);
//sL.write();

volScalarField sL0("sL0",oneMS * fH * 0.0);  // field of the correct size, filled with zeros

#include "compute_sL.H"

//dimensionedScalar sLref = 0.5*(max(sL0)+min(sL0));	// sL can become very big if Tu is high, better use sL0


volScalarField sT  // turbulent flame speed
(
    IOobject
    (
        "sT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Xi*sL,
    zeroGradientFvPatchScalarField::typeName // use zeroGradient instead of calculated BC 
);
//sT.write();



volScalarField muUnburned  // dynamic viscosity of unburned mixture
(
    IOobject
    (
        "muUnburned",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muUnburned", dimMass/dimLength/dimTime, 0.0)
);


volScalarField nuUnburned  // kinematic viscosity of unburned mixture
(
    IOobject
    (
        "nuUnburned",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("nuUnburned", dimArea/dimTime, 0.0)
);

volScalarField nu  // kinematic viscosity of mixture
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("nu", dimArea/dimTime, 0.0)
);


volScalarField alphaUnburned  // lambda/cp
(
    IOobject
    (
        "alphaUnburned",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaUnburned", dimMass/dimLength/dimTime, 0.0)
);

volScalarField aUnburned  // thermal diffusivity of unburned mixture
(
    IOobject
    (
        "aUnburned",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("aUnburned", dimArea/dimTime, 0.0)
);



volScalarField omegaC  // source term for c
(
    IOobject
    (
        "omegaC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
//    dimensionedScalar("omegaC", dimless/dimTime, 0.0)
    dimensionedScalar("omegaC", dimMass/dimVolume/dimTime, 0.0)
);
//omegaC.write();

volScalarField omegaC_ign  // source term for c
(
    IOobject
    (
        "omegaC_ign",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaC", dimMass/dimVolume/dimTime, 0.0)
);
omegaC_ign.write();

volScalarField omegaC_burn  // source term for c
(
    IOobject
    (
        "omegaC_burn",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaC", dimMass/dimVolume/dimTime, 0.0)

);
//omegaC_burn.write();

volScalarField omegaC_unquenched(omegaC_burn);


        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

	
// need inward and outward pointing face fluxes for acoustic Courant number 
surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

/*
volScalarField invDa // Damkoehler number
(
    IOobject
    (
        "invDa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("invDa", dimless, 1e+10)
);
volScalarField Da = scalar(1.0)/invDa;
//invDa.write();


volScalarField Ka // Karlovitz number
(
    IOobject
    (
        "Ka",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Ka", dimless, 1e-10)
);
//Ka.write();
*/


volScalarField Teff
(
    IOobject
    (
        "Teff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, // to continue from a stopped calculation
        IOobject::NO_WRITE
    ),
    T
);
Teff.write();


volScalarField gCr  // critical velocity gradient for quenching
(
    IOobject
    (
        "gCr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gCr", dimensionSet(0,0,-1,0,0,0,0), 1.0e5),
    zeroGradientFvPatchScalarField::typeName // use zeroGradient instead of calculated BC 
);


volScalarField G  // stretch factor for quenching
(
    IOobject
    (
        "G",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("G", dimensionSet(0,0,0,0,0,0,0), 1.0)
);
//G.write();
volScalarField epsilonCr = 0.0*turbulence->epsilon(); // needed to compute G


#include "XiEqn.H"

sL.write();
//sT.write();
omegaC_burn = uthermo.rhou() * sT * mag(fvc::grad(c)); // hier noch ohne quenching
omegaC_burn.write();
G.write();




label cIndex(0);		// these indices are used for retrieving the fluxes of the scalar quantities
label tauIndex(1);
label fHIndex(2);

UPtrList<volScalarField> rhoScalars(3);
rhoScalars.set(cIndex,&rhoC);	// this works only with a UPtrList; with a PtrList, the destructor deletes the fields it is pointing to
rhoScalars.set(tauIndex,&rhoTau);
rhoScalars.set(fHIndex,&rhofH);



    Info<< "\nCreate Riemann solver: "; // Roe does not work yet
    
     RiemannFlux<hllcALEFlux> Riemann(rhoU, rho, rhoE, rhoEu, rhoScalars, thermo, turbulence());	word RiemannSolver("HLLC");
     // for the other schemes, rhoEu still needs to be implemented
     //RiemannFlux<AUSMplusALEFlux> Riemann(rhoU, rho, rhoE, rhoScalars, thermo, turbulence());	word RiemannSolver("AUSM+");
     //RiemannFlux<AUSMplusUpALEFlux> Riemann(rhoU, rho, rhoE, rhoScalars, thermo, turbulence());	word RiemannSolver("AUSM+Up");
     //RiemannFlux<rusanovALEFlux>  Riemann(rhoU, rho, rhoE, rhoScalars, thermo, turbulence());	word RiemannSolver("Rusanov");

     Info << RiemannSolver << endl;
     


volScalarField cSound
(
    IOobject
    (
        "cSound",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt((thermo.Cp()/thermo.Cv())/thermo.psi()),
    zeroGradientFvPatchScalarField::typeName // use zeroGradient instead of calculated BC , better for post-processing
);
cSound.write();


// MinMax model:
dimensionedScalar dTmin("dTmin",dimTemperature,scalar(0.1));
dimensionedScalar dpmin("dpmin",dimPressure,scalar(10.0));

volScalarField alphaHI
(
    IOobject
    (
        "alphaHI",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);

volScalarField pHI
(
    IOobject
    (
        "pHI",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    p+dpmin
);


volScalarField pLO
(
    IOobject
    (
        "pLO",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    p-dpmin
);


volScalarField THI
(
    IOobject
    (
        "THI",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    T+dTmin
);


volScalarField TLO
(
    IOobject
    (
        "TLO",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    T-dTmin
);

volScalarField tauHI
(
    IOobject
    (
        "tauHI",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    tau
);
tauHI.max(0.0);
tauHI.min(1.0);
//tauHI.write();

volScalarField tauLO
(
    IOobject
    (
        "tauLO",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    tau
);
tauLO.max(0.0);
tauLO.min(1.0);
//tauLO.write();

volScalarField omegaTauHI(omegaTau);
volScalarField omegaTauLO(omegaTau);
volScalarField tIgnHI(tIgn);
volScalarField tIgnLO(tIgn);



if(reconstructionStages)
{
#include "findMinMax_fromCellAvg.H"
}
/*
alphaHI.write();
pHI.write();
pLO.write();
TeffHI.write();
TeffLO.write();
*/

p.oldTime();

Info << "createFields.H complete" << endl;
