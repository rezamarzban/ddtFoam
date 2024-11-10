Info<< "Reading additional properties\n" << endl;

IOdictionary additionalProperties
(
    IOobject
    (
        "additionalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

Switch inviscid(additionalProperties.lookup("inviscid"));
Info << "\n --- The simulation is ";
	if (inviscid) Info << "inviscid";
	else Info << "viscous !!!";
	Info << " ---\n";
Info << nl << endl;


scalar ScT(additionalProperties.lookupOrDefault("ScT",1.0)); // turbulent Schmidt number
Info << "The turbulent Schmidt number is " << ScT << endl;

Switch temperatureFix(additionalProperties.lookupOrDefault("temperatureFix",false));
//Info << "The temperature and density fix is " << temperatureFix << endl;


