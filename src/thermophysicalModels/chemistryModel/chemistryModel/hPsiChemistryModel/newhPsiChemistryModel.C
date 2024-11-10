/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "hPsiChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::hPsiChemistryModel> Foam::hPsiChemistryModel::New
(
    const fvMesh& mesh
)
{
    word hPsiChemistryModelType;
    word thermoTypeName;
    word userModel;

    // Enclose the creation of the chemistrtyProperties to ensure it is
    // deleted before the chemistrtyProperties is created otherwise the
    // dictionary is entered in the database twice
    {
        IOdictionary chemistryPropertiesDict
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        chemistryPropertiesDict.lookup("hPsiChemistryModel") >> userModel;

        // construct chemistry model type name by inserting first template
        // argument
        label tempOpen = userModel.find('<');
        label tempClose = userModel.find('>');

        word className = userModel(0, tempOpen);
        thermoTypeName = userModel(tempOpen + 1, tempClose - tempOpen - 1);

        hPsiChemistryModelType =
            className + '<' + typeName + ',' + thermoTypeName + '>';
    }

    if (debug)
    {
        Info<< "Selecting hPsiChemistryModel " << hPsiChemistryModelType << endl;
    }
    else
    {
        Info<< "Selecting hPsiChemistryModel " << userModel << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(hPsiChemistryModelType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        if (debug)
        {
            FatalErrorIn("hPsiChemistryModelBase::New(const mesh&)")
                << "Unknown hPsiChemistryModel type " << hPsiChemistryModelType
                << nl << nl << "Valid hPsiChemistryModel types are:" << nl
                << fvMeshConstructorTablePtr_->toc() << nl << exit(FatalError);
        }
        else
        {
            wordList models = fvMeshConstructorTablePtr_->toc();
            forAll(models, i)
            {
                models[i] = models[i].replace(typeName + ',', "");
            }

            FatalErrorIn("hPsiChemistryModelBase::New(const mesh&)")
                << "Unknown hPsiChemistryModel type " << userModel
                << nl << nl << "Valid hPsiChemistryModel types are:" << nl
                << models << nl << exit(FatalError);
        }
    }

    return autoPtr<hPsiChemistryModel>
        (cstrIter()(mesh, typeName, thermoTypeName));
}


// ************************************************************************* //
