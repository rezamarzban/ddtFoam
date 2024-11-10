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

#include "gradientUnburnedEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "unburnedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientUnburnedEnthalpyFvPatchScalarField::gradientUnburnedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::gradientUnburnedEnthalpyFvPatchScalarField::gradientUnburnedEnthalpyFvPatchScalarField
(
    const gradientUnburnedEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::gradientUnburnedEnthalpyFvPatchScalarField::gradientUnburnedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::gradientUnburnedEnthalpyFvPatchScalarField::gradientUnburnedEnthalpyFvPatchScalarField
(
    const gradientUnburnedEnthalpyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::gradientUnburnedEnthalpyFvPatchScalarField::gradientUnburnedEnthalpyFvPatchScalarField
(
    const gradientUnburnedEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientUnburnedEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

      gradient() = uthermo.Cpu(Tuw, fHw, patchi)*Tuw.snGrad()
        + patch().deltaCoeffs()*
        (
            uthermo.hu(Tuw, fHw, patchi)
          - uthermo.hu(Tuw, fHw, patch().faceCells())
        );
  
/*
    if (dimensionedInternalField().name() == "h")
    {
        gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
        + patch().deltaCoeffs()*
        (
            thermo.h(Tw, patchi)
          - thermo.h(Tw, patch().faceCells())
        );
    }
    else
    {
        gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
        + patch().deltaCoeffs()*
        (
            thermo.hs(Tw, patchi)
          - thermo.hs(Tw, patch().faceCells())
        );
    }
*/
    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradientUnburnedEnthalpyFvPatchScalarField
    );
}


// ************************************************************************* //
