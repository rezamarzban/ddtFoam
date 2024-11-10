const volScalarField::GeometricBoundaryField& tbf = T.boundaryField();
wordList rhoEBoundaryTypes = tbf.types();

forAll(rhoEBoundaryTypes, patchi)
{

      if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
      {
          rhoEBoundaryTypes[patchi] = fixedRhoEFvPatchScalarField::typeName;

      }
	// still to be implemented:
      /*
      else if(isA<zeroGradientFvPatchScalarField>(tbf[patchi]) || isA<fixedGradientFvPatchScalarField>(tbf[patchi]))
      {
	  rhoEBoundaryTypes[patchi] = gradientrhoEFvPatchScalarField::typeName;
      }
      else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
      {
	  rhoEBoundaryTypes[patchi] = mixedrhoEFvPatchScalarField::typeName;
      }
      */
      else //if(isA<zeroGradientFvPatchScalarField>(tbf[patchi]))
      {
          rhoEBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;

      }
}
