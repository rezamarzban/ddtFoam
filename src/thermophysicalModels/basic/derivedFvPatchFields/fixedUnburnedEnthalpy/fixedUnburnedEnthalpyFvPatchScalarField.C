/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "fixedUnburnedEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "unburnedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedUnburnedEnthalpyFvPatchScalarField::fixedUnburnedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::fixedUnburnedEnthalpyFvPatchScalarField::fixedUnburnedEnthalpyFvPatchScalarField
(
    const fixedUnburnedEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixedUnburnedEnthalpyFvPatchScalarField::fixedUnburnedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::fixedUnburnedEnthalpyFvPatchScalarField::fixedUnburnedEnthalpyFvPatchScalarField
(
    const fixedUnburnedEnthalpyFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::fixedUnburnedEnthalpyFvPatchScalarField::fixedUnburnedEnthalpyFvPatchScalarField
(
    const fixedUnburnedEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedUnburnedEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//	Info << db().names() << endl;

    const unburnedThermo& uthermo = db().lookupObject<unburnedThermo>
    (
        "uthermoProperties"
    );

    const label patchi = patch().index();

    fvPatchScalarField& Tuw =
        const_cast<fvPatchScalarField&>(uthermo.Tu().boundaryField()[patchi]);
    Tuw.evaluate();

    fvPatchScalarField& fHw =
        const_cast<fvPatchScalarField&>(uthermo.fH().boundaryField()[patchi]);
    fHw.evaluate();


        operator==(uthermo.hu(Tuw, fHw, patchi));

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedUnburnedEnthalpyFvPatchScalarField
    );
}


// ************************************************************************* //
