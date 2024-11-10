
// find min and max values from averaged values of the surrounding cells


volScalarField kappa = thermo.Cp() / thermo.Cv();

forAll(mesh.cells(),celli)
{
      
	
	DynamicList<scalar> pNList;
	pNList.append(p[celli]);
	
	if(reconstructionStages>0)
	{
	    labelList neighbourList = mesh.cellCells()[celli];
	
	    forAll(neighbourList,nei)
	    {
		pNList.append(p[neighbourList[nei]]);
	    }
	}
   
    pLO[celli] = max(0.9e5,min(pNList));	// bounding
    pHI[celli] = max(max(pNList),pLO[celli]+dpmin.value());
    
    alphaHI[celli] = (p[celli]-pLO[celli])/(pHI[celli]-pLO[celli]);
}


//Info<< "  kappa   min...max   = " << min(kappa).value() << " ... " << max(kappa).value() << endl;
volScalarField MaShock = sqrt((pHI/pLO - scalar(1.0))*(kappa+scalar(1.0))/(2*kappa)+scalar(1.0));
//Info<< "  MaShock   min...max   = " << min(MaShock).value() << " ... " << max(MaShock).value() << endl;
volScalarField Tratio = pHI/pLO * (scalar(1.0)-scalar(2.0)/(kappa+scalar(1.0))*(scalar(1.0)-pow(MaShock,-2.0)));
//Info<< "  Tratio   min...max   = " << min(Tratio).value() << " ... " << max(Tratio).value() << endl;

TLO = T / (alphaHI*Tratio+scalar(1.0)-alphaHI);
THI = Tratio * TLO;



    
    
