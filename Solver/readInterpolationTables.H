Info<< "Reading interpolation tables" << endl;

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

IOdictionary tignDictionary
(
    IOobject
    (
        "tignTable",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);



interpolationLookUpTable<scalar> cTable(cDictionary);
interpolationLookUpTable<scalar> tignTable(tignDictionary);

Info << "tignTable is logarithmic: " << tignTable.islogarithmic() << endl;

if(!tignTable.islogarithmic())
{
       FatalError << args.executable() << " : provide a logarithmic tignTable! " << abort(FatalError);   
}


