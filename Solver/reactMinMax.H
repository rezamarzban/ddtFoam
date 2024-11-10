{
    omegaC_ign *= 0.0;
 
    
    if(autoIgnitionTauModel)
    { 
      
      omegaC_ign += alphaHI * rho * (scalar(1.0)-c) / runTime.deltaT() * Foam::pos(tauHI-scalar(0.99));
      omegaC_ign += (scalar(1.0)-alphaHI) * rho * (scalar(1.0)-c) / runTime.deltaT() * Foam::pos(tauLO-scalar(0.99));	

      omegaC_ign = min(omegaC_ign,rho*(scalar(1.0)-c)/runTime.deltaT());
      omegaC_ign.max(0.0);
      
    }
    
    
    const label tautableDim = 3; 
    
    forAll(p,i)
    {
	List<scalar> lookupCoordinates(tautableDim,0.0);
		lookupCoordinates[0]=fH[i];
		lookupCoordinates[1]=pHI[i];		// old values, no thermo.update() 
		lookupCoordinates[2]=THI[i];

    	const List<scalar> &LOG_t_ign_temp = tignTable.lookUpBounded(lookupCoordinates,false);
	tIgnHI[i] = min(max(1e-14,pow(10,LOG_t_ign_temp[0])),1e9);
	
    }
    tIgnHI.correctBoundaryConditions(); // only for post-processing: update zeroGradient BC
    
    omegaTauHI = rho / tIgnHI;
    
    forAll(p,i)
    {
	List<scalar> lookupCoordinates(tautableDim,0.0);
		lookupCoordinates[0]=fH[i];
		lookupCoordinates[1]=pLO[i];		// old values, no thermo.update() 
		lookupCoordinates[2]=TLO[i];

    	const List<scalar> &LOG_t_ign_temp = tignTable.lookUpBounded(lookupCoordinates,false);
	tIgnLO[i] = min(max(1e-14,pow(10,LOG_t_ign_temp[0])),1e9);
    }
    tIgnLO.correctBoundaryConditions(); // only for post-processing: update zeroGradient BC
    
    omegaTauLO = rho / tIgnLO;
    
    
    omegaTau = alphaHI * omegaTauHI + (scalar(1.0)-alphaHI) * omegaTauLO;
    
    tIgn = tIgnHI; // for post-processing
        
}
