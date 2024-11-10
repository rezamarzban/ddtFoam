
 
const dimensionedScalar oneMS("oneMS", dimensionSet(0,1,-1,0,0), 1.0); // for sL correlation
const dimensionedScalar oneM("oneM", dimensionSet(0,1,0,0,0), 1.0); // for Zbikowski correlation
const dimensionedScalar oneKGM3("oneKGM3", dimensionSet(1,-3,0,0,0), 1.0); // for rhoUnburned
const dimensionedScalar oneM2S("oneM2S", dimensionSet(0,2,-1,0,0), 1.0); // for aUnburned / nuUnburned
const dimensionedScalar oneKGMS("oneKGMS", dimensionSet(1,-1,-1,0,0), 1.0); // for alphaUnburned / etaUnburned
const dimensionedScalar pref("pref", dimensionSet(1,-1,-2,0,0), 1.0e5); // for rhoUnburned/alphaUnburned/etaUnburned
const dimensionedScalar pmin("pmin", dimensionSet(1,-1,-2,0,0), 0.2e5); // for bounding


Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<hPsiChemistryModel> pChemistry
(
    hPsiChemistryModel::New(mesh)
);


hPsiChemistryModel& chemistry = pChemistry();

hCombustionThermo& thermo = chemistry.thermo();


basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

Info << "thermo: " << thermo << endl; 
Info << "species: " << composition.species() << endl; 


word inertSpecie(thermo.lookup("inertSpecie"));

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
    thermo.rho()
);
rho.write();

Info<< "Reading field U\n" << endl;
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


volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& h = thermo.h();
//volScalarField& e = thermo.e();

const volScalarField& T = thermo.T();



//Info<< "h  min...max   = " << min(h).value() << " ... " << max(h).value() << endl;

#include "compressibleCreatePhi.H"


volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

Info << "Creating turbulence model.\n" << nl;
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

Info<< "Creating field DpDt\n" << endl;
//volScalarField DpDt =    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

volScalarField DpDt
(
    IOobject
    (
        "DpDt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
);



multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

label h2Index=-1;
label h2oIndex=-1;
label o2Index=-1;
label n2Index=-1;
forAll(Y, i)
{
    fields.add(Y[i]);
	if(Y[i].name()=="H2") h2Index=i;
	if(Y[i].name()=="H2O") h2oIndex=i;
	if(Y[i].name()=="O2") o2Index=i;
	if(Y[i].name()=="N2") n2Index=i;
}
// alternativ:
// h2Index = composition.species()["H2"];


/*
volScalarField& yH2  = Y[h2Index];
volScalarField& yH2O = Y[h2oIndex];
volScalarField& yO2  = Y[o2Index];
volScalarField& yN2  = Y[n2Index];
*/


// Info << "Indices: \nH2:  " << h2Index <<"\nO2:  " << o2Index << "\nN2:  " << n2Index << "\nH2O: " << h2oIndex << nl << endl;





//fields.add(e);
fields.add(h);



volScalarField Cp = thermo.Cp();
volScalarField Cv = thermo.Cv();
// first had to implement Cv in
// src/thermophysicalModels/reactionThermo/combustionThermo/mixtureThermos/hsPsiMixtureThermo 


volScalarField cSound  // speed of sound
(
    IOobject
    (
        "cSound",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
//    mesh,
//    dimensionedScalar("zero", dimensionSet(0,1,-1,0,0,0,0), 0.0)
    sqrt((Cp/Cv)*(Cp - Cv)*T)
);
cSound.write();	


volScalarField c  // reaction progress
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

#include "createPhiC.H"




volScalarField fH  // H2 mixture fraction (unburned)
(
    IOobject
    (
        "fH",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Y[h2Index] + MH2/MH2O * Y[h2oIndex],
    Y[h2Index].boundaryField().types()
);
//fH.boundaryField() = Y[h2Index].boundaryField() + MH2/MH2O * Y[h2oIndex].boundaryField();
fH.write();


    unburnedThermo uthermo(mesh, p, fH, yO2inAir);
    const volScalarField& Tu=uthermo.Tu();
    Info << " Tu varies between " << min(Tu).value() << " and " << max(Tu).value() << " K." << endl;
    const volScalarField& rhou = uthermo.rhou();
    volScalarField& hu = uthermo.hu();		


volScalarField omegaC  // source term for c
(
    IOobject
    (
        "omegaC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaC", dimMass/dimVolume/dimTime, 0.0)
);
omegaC.write();


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
    dimensionedScalar("sL", dimLength/dimTime, 0.0)
);
//sL.write();


volScalarField sT  // turbulent flame speed
(
    IOobject
    (
        "sT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sL
);
//sT.write();



volScalarField rhoUnburned  // density of unburned mixture
(
    IOobject
    (
        "rhoUnburned",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoUnburned", dimMass/dimVolume, 0.0)
);


volScalarField etaUnburned  // dynamic viscosity of unburned mixture
(
    IOobject
    (
        "etaUnburned",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("etaUnburned", dimMass/dimLength/dimTime, 0.0)
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

volScalarField epsilonCr = 0.0*turbulence->epsilon();

volScalarField G  // stretch factor
(
    IOobject
    (
        "G",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("G", dimensionSet(0,0,0,0,0,0,0), 1.0)
);
//G.write();

volScalarField DpDtold = DpDt;


volScalarField Da // Damkoehler number
(
    IOobject
    (
        "Da",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Da", dimless, 0.0)
);
//Da = pow(sL,2)/nuUnburned*turbulence->k()/turbulence->epsilon();
//Da.write();

/*
thermo.correct(); // update es, h and hs for postprocessing
*/

/*

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );


volScalarField gradC
(
    IOobject
    (
        "gradC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gradC", dimensionSet(0,-1,0,0,0), 0.0)   
);

volScalarField Re_t	// turbulent Reynolds number, updated in tfcSource.H
(
    IOobject
    (
        "Re_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Re_t", dimless, 0.0)   
);
*/
/*
volScalarField UpBySl	// u prime divided by laminar flame speed
(
    IOobject
    (
        "UpBySl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("UpBySl", dimless, 0.0)   
);
*/




volScalarField sL0("sL0",oneMS * fH * 0.0);  // field of the correct size, filled with zeros

#include "compute_sL.H"

dimensionedScalar sLref = 0.5*(max(sL)+min(sL));


    volScalarField b
    (
        IOobject
        (
            "b",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        c // preliminary: takes boundary conditions from c
    );
    b = scalar(1) - c;
 

    Info<< "Creating field Xi\n" << endl;
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




Info << endl;



#include "acousticCreatePhi.H"


