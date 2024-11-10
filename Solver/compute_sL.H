{

if(sLCorrelation=="Ardey")	// exp. data from Koroll et al., 1993; stretch effects have not been considered
{
    sL0 = oneMS * (-40452.0 * pow(fH,4) + 21026.0 * pow(fH,3) - 3721.6 * pow(fH,2) + 238.11 * fH -1.5938);
}
else //if(sLCorrelation=="Ettner")	// exp. data from Konnov, 2008; considered stretch effects in experiments
{
   // unburned molar H fraction:
   volScalarField xH("xH", (fH/MH2) / (fH/MH2 + yO2inAir*(scalar(1.0)-fH)/MO2 + (1.0-yO2inAir)*(scalar(1.0)-fH)/MN2)); 
   
   sL0 = oneMS * 
	 (
	    (-4.8889e+02 * pow(xH,4) + 2.8497e+02 * pow(xH,3) - 2.1924e+01 * pow(xH,2) + 1.3506e+00 * xH -4.0000e-02) *Foam::pos(0.35-xH) 
	   +(-1.6016e+02 * pow(xH,4) + 3.7773e+02 * pow(xH,3) - 3.4865e+02 * pow(xH,2) + 1.3996e+02 * xH -1.7448e+01) *Foam::pos(xH-0.35) 
	 );
}

  const dimensionedScalar alpha("alpha", dimless, 1.75);
  const dimensionedScalar beta("beta", dimless, -0.2);
  const dimensionedScalar Tref("Tref", dimensionSet(0,0,0,1,0), 298.0);
  const dimensionedScalar pRef("pRef", dimensionSet(1,-1,-2,0,0), 1.013e5);
  sL = sL0 * pow((Tu/Tref), alpha) * pow((p/pRef), beta);
  sL.correctBoundaryConditions(); // to avoid 0 on BCs, required for divisions by sL
  sL.max(sLMin);
}
