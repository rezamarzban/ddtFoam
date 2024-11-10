{
//Info << "looking up y from c" << endl;

const label ctableDim = 3; 


forAll(c,i)
{
  
	const double& fi = fH[i];

	List<scalar> yunburned(Y.size(),0.0);
	    yunburned[h2Index]=fH[i];
	    yunburned[o2Index]=yO2inAir * (1.0-fi);;
	    yunburned[n2Index]=1.0-sum(yunburned);


	List<scalar> lookupCoordinates(ctableDim,0.0);
	
	
		lookupCoordinates[0]=fH[i];
		lookupCoordinates[1]=p[i];	
		lookupCoordinates[2]=Tu[i];

	const List<scalar> &yburned = cTable.lookUp(lookupCoordinates,false);

	forAll(Y,j)
	{
	    Y[j][i] = (1.0-c[i])*yunburned[j] + c[i]*yburned[j]; 	
	}


}

forAll(Y,j)
{
	Y[j].correctBoundaryConditions();	// required because no transport eq. is solved for the species
}	

}