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

#include "mixedUnburnedEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "eUnburnedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedUnburnedEnergyFvPatchScalarField::mixedUnburnedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixedUnburnedEnergyFvPatchScalarField::mixedUnburnedEnergyFvPatchScalarField
(
    const mixedUnburnedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedUnburnedEnergyFvPatchScalarField::mixedUnburnedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedUnburnedEnergyFvPatchScalarField::mixedUnburnedEnergyFvPatchScalarField
(
    const mixedUnburnedEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixedUnburnedEnergyFvPatchScalarField::mixedUnburnedEnergyFvPatchScalarField
(
    const mixedUnburnedEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedUnburnedEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const eUnburnedThermo& uthermo = db().lookupObject<eUnburnedThermo>
    (
        "uthermoProperties"
    );

    const label patchi = patch().index();

    mixedFvPatchScalarField& Tuw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(uthermo.Tu().boundaryField()[patchi])
    );
    Tuw.evaluate();

    mixedFvPatchScalarField& fHw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(uthermo.fH().boundaryField()[patchi])
    );
    fHw.evaluate();



    valueFraction() = Tuw.valueFraction();

        refValue() = uthermo.eu(Tuw.refValue(), fHw.refValue(), patchi);
        refGrad() = uthermo.Cvu(Tuw, fHw, patchi)*Tuw.refGrad()
        + patch().deltaCoeffs()*
         (
            uthermo.eu(Tuw, fHw, patchi)
          - uthermo.eu(Tuw, fHw, patch().faceCells())
         );



/*
    if (dimensionedInternalField().name() == "h")
    {
        refValue() = thermo.h(Tw.refValue(), patchi);
        refGrad() = thermo.Cp(Tw, patchi)*Tw.refGrad()
        + patch().deltaCoeffs()*
         (
            thermo.h(Tw, patchi)
          - thermo.h(Tw, patch().faceCells())
         );
    }
    else
    {
        refValue() = thermo.hs(Tw.refValue(), patchi);
        refGrad() = thermo.Cp(Tw, patchi)*Tw.refGrad()
        + patch().deltaCoeffs()*
         (
            thermo.hs(Tw, patchi)
          - thermo.hs(Tw, patch().faceCells())
         );
    }
*/
    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedUnburnedEnergyFvPatchScalarField
    );
}


// ************************************************************************* //
