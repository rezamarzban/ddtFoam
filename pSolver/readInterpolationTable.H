Info<< "Reading interpolation table for mixture composition" << endl;

IOdictionary cDictionary
(
    IOobject
    (
        "cTable",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

interpolationLookUpTable<scalar> cTable(cDictionary);
