const volScalarField::GeometricBoundaryField& pbf = p.boundaryField();
const volScalarField::GeometricBoundaryField& Tbf = T.boundaryField();
wordList rhoBoundaryTypes = pbf.types();	// take the same BCs as pressure

forAll(rhoBoundaryTypes, patchi)	// except:
{
    if (rhoBoundaryTypes[patchi] == "waveTransmissive")
    {
        rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
    else if (pbf[patchi].fixesValue())
    {
        rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField::typeName;
    }    
    else if (Tbf[patchi].fixesValue())
    {
        rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField::typeName;
    }
/*    
    else
    {
        rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
*/
}
