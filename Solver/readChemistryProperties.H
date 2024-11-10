
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

scalar yO2inAir(0.23295);
scalar reconstructionStages(1);




    word sLCorrelation("Ettner");    
    word XiModel("transport");
    dimensionedScalar XiCoef("XiCoef", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.62);
    dimensionedScalar XiShapeCoef("XiShapeCoef", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1);
    dimensionedScalar uPrimeCoef("uPrimeCoef", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1);
    
Switch turbulentBurning(chemistryProperties.lookup("turbulentBurning"));
Switch autoIgnitionTauModel(chemistryProperties.lookup("autoIgnitionTauModel"));

    reconstructionStages = readScalar(chemistryProperties.lookup("reconstructionStages"));
    Info << "  reconstructionStages = " << reconstructionStages << endl;


dictionary tfcDict = chemistryProperties.subDict("tfcDict");
    yO2inAir = readScalar(tfcDict.lookup("yO2inAir"));
    
	sLCorrelation = word(tfcDict.lookupOrAddDefault<word>("sLCorrelation","Ettner",false,true));
        XiModel = word(tfcDict.lookup("XiModel"));
        XiCoef = tfcDict.lookup("XiCoef");
        XiShapeCoef = tfcDict.lookup("XiShapeCoef");
        uPrimeCoef = tfcDict.lookup("uPrimeCoef");
        Info << "  TFC model: bXi with" << endl;
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
	
	Switch quenching(chemistryProperties.lookupOrDefault("quenching",true));
	Info << "\nQuenching is " << quenching << endl;
