Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);



dictionary tfcDict = chemistryProperties.subDict("tfcDict");

    scalar yO2inAir(0.23295);

    // heat of formation of H2O; J/kmol:
    const dimensionedScalar hf_H2O("hf_H2O", dimensionSet(1,2,-2,0,-1),-241845e3); 

    word speciesSolver("transport");

    word tfcModel("Zimont");

    Switch quenching(chemistryProperties.lookup("quenching"));

    Info << "Quenching is " ;
        if(!quenching) Info << "NOT ";
    Info << "taken into account" << endl;

    word sLCorrelation("Ettner");
    word XiModel("transport");
    dimensionedScalar XiCoef("XiCoef", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.62);
    dimensionedScalar XiShapeCoef("XiShapeCoef", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1);
    dimensionedScalar uPrimeCoef("uPrimeCoef", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1);


    Info << nl << " -- Using the b-Xi model -- " << endl;

    tfcModel = "bXi";
	sLCorrelation = word(tfcDict.lookupOrAddDefault<word>("sLCorrelation","Ettner",false,true));
        XiModel = word(tfcDict.lookup("XiModel"));
        XiCoef = tfcDict.lookup("XiCoef");
        XiShapeCoef = tfcDict.lookup("XiShapeCoef");
        uPrimeCoef = tfcDict.lookup("uPrimeCoef");
        Info << "TFC model: bXi with" << endl;
	if((sLCorrelation=="Ettner") || (sLCorrelation=="Ardey"))
	  { Info << "  sLCorrelation = " << sLCorrelation << endl;}
	else
	  {       FatalError
		  << args.executable() << " : unknown sLCorrelation " << sLCorrelation
		  << abort(FatalError);  
	  }
        Info << "  XiModel = " << XiModel << endl;
        Info << "  XiCoef = " << XiCoef.value() << endl;
        Info << "  XiShapeCoef = " << XiShapeCoef.value() << endl;
        Info << "  uPrimeCoef = " << uPrimeCoef.value() << endl;



    yO2inAir = readScalar(tfcDict.lookup("yO2inAir"));
    Info << "mass fraction of O2 in air = " << yO2inAir << endl;


    Info << endl;


Info << endl;