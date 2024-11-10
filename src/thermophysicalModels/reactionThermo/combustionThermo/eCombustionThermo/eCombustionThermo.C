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

#include "eCombustionThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eCombustionThermo, 0);
    defineRunTimeSelectionTable(eCombustionThermo, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eCombustionThermo::eCombustionThermo(const fvMesh& mesh)
:
    basicPsiThermo(mesh),

    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->eBoundaryTypes()
    ),

    es_
    (
        IOobject
        (
            "es",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->eBoundaryTypes()
    ),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->hBoundaryTypes()
    ),
    hs_
    (
        IOobject
        (
            "hs",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->hBoundaryTypes()
    )


{
	Info << "eCombustionThermo created" << endl;
	//Info << " (values of e_@fixedValue are still zero)" << endl;
	//this->eBoundaryCorrection(e_);	// bringt nix, bleiben null
	//Info << "corrected e_: " << e_.boundaryField() << endl;
}
/*
Foam::tmp<Foam::volScalarField> Foam::eCombustionThermo::e
(
    volScalarField& Tgiven
) const
{
    notImplemented
    (
        "basicThermo::e"
        "(volScalarField& Tgiven) const"
    );
    return tmp<volScalarField>(NULL);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eCombustionThermo::~eCombustionThermo()
{}


// ************************************************************************* //
/*
Foam::volScalarField Foam::eCombustionThermo::eFromT(const volScalarField& Tgiven) const
{
  //Info << " eCombustionThermo::eFromT " << endl;
    notImplemented("eCombustionThermo::eFromT(volScalarField& Tgiven) const");
    return volScalarField::null();
  //return e_;
}
*/
Foam::tmp<Foam::scalarField> Foam::eCombustionThermo::e
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "eCombustionThermo::e"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(NULL);
}


void Foam::eCombustionThermo::eFromT(volScalarField& eToCompute, const volScalarField& Tgiven) const
{
    notImplemented("eCombustionThermo::eFromT(volScalarField& eToCompute, const volScalarField& Tgiven) const");
}
void Foam::eCombustionThermo::hFromT(volScalarField& hToCompute, const volScalarField& Tgiven) const
{
    notImplemented("eCombustionThermo::hFromT(volScalarField& hToCompute, const volScalarField& Tgiven) const");
}

