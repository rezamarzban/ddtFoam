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

#include "mixedUnburnedEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "unburnedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedUnburnedEnthalpyFvPatchScalarField::mixedUnburnedEnthalpyFvPatchScalarField
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


Foam::mixedUnburnedEnthalpyFvPatchScalarField::mixedUnburnedEnthalpyFvPatchScalarField
(
    const mixedUnburnedEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedUnburnedEnthalpyFvPatchScalarField::mixedUnburnedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedUnburnedEnthalpyFvPatchScalarField::mixedUnburnedEnthalpyFvPatchScalarField
(
    const mixedUnburnedEnthalpyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixedUnburnedEnthalpyFvPatchScalarField::mixedUnburnedEnthalpyFvPatchScalarField
(
    const mixedUnburnedEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedUnburnedEnthalpyFvPatchScalarField::updateCoeffs()
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

        refValue() = uthermo.hu(Tuw.refValue(), fHw.refValue(), patchi);
        refGrad() = uthermo.Cpu(Tuw, fHw, patchi)*Tuw.refGrad()
        + patch().deltaCoeffs()*
         (
            uthermo.hu(Tuw, fHw, patchi)
          - uthermo.hu(Tuw, fHw, patch().faceCells())
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
        mixedUnburnedEnthalpyFvPatchScalarField
    );
}


// ************************************************************************* //
