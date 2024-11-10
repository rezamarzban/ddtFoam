    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() == inertSpecie) inertIndex = i;
	else
	{
	    Yt += Y[i];
	}
    }
    Info << "no transport equation solved for specie: Y[" << inertIndex << "] = " << Y[inertIndex].name() << endl;
     
    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

