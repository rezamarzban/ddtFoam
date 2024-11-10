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

Foam::autoPtr<Foam::eCombustionThermo> Foam::eCombustionThermo::New
(
    const fvMesh& mesh
)
{
    word eCombustionThermoTypeName;

    // Enclose the creation of the thermophysicalProperties to ensure it is
    // deleted before the turbulenceModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        thermoDict.lookup("thermoType") >> eCombustionThermoTypeName;
    }

    Info<< "Selecting thermodynamics package " << eCombustionThermoTypeName
        << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(eCombustionThermoTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("eCombustionThermo::New(const fvMesh&)")
            << "Unknown eCombustionThermo type "
            << eCombustionThermoTypeName << nl << nl
            << "Valid eCombustionThermo types are:" << nl
            << fvMeshConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<eCombustionThermo>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::eCombustionThermo> Foam::eCombustionThermo::NewType
(
    const fvMesh& mesh,
    const word& thermoType
)
{
    word eCombustionThermoTypeName;

    // Enclose the creation of the thermophysicalProperties to ensure it is
    // deleted before the turbulenceModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        thermoDict.lookup("thermoType") >> eCombustionThermoTypeName;

        if (eCombustionThermoTypeName.find(thermoType) == string::npos)
        {
            wordList allModels = fvMeshConstructorTablePtr_->sortedToc();
            DynamicList<word> validModels;
            forAll(allModels, i)
            {
                if (allModels[i].find(thermoType) != string::npos)
                {
                    validModels.append(allModels[i]);
                }
            }

            FatalErrorIn
            (
                "autoPtr<eCombustionThermo> eCombustionThermo::NewType"
                "("
                    "const fvMesh&, "
                    "const word&"
                ")"
            )   << "Inconsistent thermo package selected:" << nl << nl
                << eCombustionThermoTypeName << nl << nl << "Please select a "
                << "thermo package based on " << thermoType
                << ". Valid options include:" << nl << validModels << nl
                << exit(FatalError);
        }
    }

    Info<< "Selecting thermodynamics package " << eCombustionThermoTypeName
        << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(eCombustionThermoTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("eCombustionThermo::New(const fvMesh&)")
            << "Unknown eCombustionThermo type "
            << eCombustionThermoTypeName << nl << nl
            << "Valid eCombustionThermo types are:" << nl
            << fvMeshConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<eCombustionThermo>(cstrIter()(mesh));
}


// ************************************************************************* //
