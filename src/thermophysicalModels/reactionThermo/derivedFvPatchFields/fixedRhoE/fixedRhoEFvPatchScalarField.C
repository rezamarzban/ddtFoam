/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedRhoEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "eCombustionThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedRhoEFvPatchScalarField::
fixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::fixedRhoEFvPatchScalarField::
fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixedRhoEFvPatchScalarField::
fixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::fixedRhoEFvPatchScalarField::
fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::fixedRhoEFvPatchScalarField::
fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedRhoEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const eCombustionThermo& thermo = db().lookupObject<eCombustionThermo>
    (
        "thermophysicalProperties"
    );

    const label patchi = patch().index();

    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    Tw.evaluate();
    
    // gets rho on boundary from fixedRhoFvPatchScalarField (rho = p * psi, psi from thermo.correct())
    // for rho to be correct, the rhoEqn needs to be solved before the rhoE Eqn.
    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>("rho");
	
    const fvPatchField<vector>& Up = 
	patch().lookupPatchField<volVectorField, vector>("U");
    
	
    	  //Info << "U          on patch" << patchi << ": " << Up << endl;
	  //Info << "U patchField on patch" << patchi << ": " << Up.Field() << endl;
	  //Info << "U patchInternalField on patch" << patchi << ": " << Up.patchInternalField() << endl;
	  //Info << "ekin of patchInternalField" << patchi << ": " << 0.5*magSqr(Up.patchInternalField()) << endl;	  
	  //Info << "ekin directly on patch = " << 0.5*magSqr(Up) << endl;
	  
    operator==(rhop * (thermo.e(Tw, patchi) + 0.5*magSqr(Up)));
    //Info << "rhoE on patch " << patchi << ": " << endl;
    //Info << rhop << nl << " < * > " << nl << thermo.e(Tw,patchi) << nl << " < + > " << nl << 0.5*magSqr(Up) << nl << " = " << nl << rhop * (thermo.e(Tw, patchi) + 0.5*magSqr(Up)) << endl;
    

    // Alternative:
    // take internal fields of rho and U - a difference of rho and U between wall and flow should not affect the rhoE eqn.
    // only a difference in e due to a temperature difference needs to be evaluated    
    /*
    //operator==(rhop.patchInternalField() * (thermo.e(Tw, patchi) + 0.5*magSqr(Up.patchInternalField())));
    operator==(rhop.patchInternalField() * (thermo.e(Tw, patchi) + 0.5*magSqr(Up)));
        //Info << "rhoE on patch " << patchi << ": " << endl;
    //Info << rhop.patchInternalField() << nl << " < * > " << nl << thermo.e(Tw,patchi) << nl << " < + > " << nl << 0.5*magSqr(Up.patchInternalField()) << nl << " = " << nl << rhop.patchInternalField() * (thermo.e(Tw, patchi) + 0.5*magSqr(Up.patchInternalField())) << endl;
    */
    
    
    

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedRhoEFvPatchScalarField
    );
}

// ************************************************************************* //
