/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "ePsiChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ePsiChemistryModel> Foam::ePsiChemistryModel::New
(
    const fvMesh& mesh
)
{
    word ePsiChemistryModelType;
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

        chemistryPropertiesDict.lookup("ePsiChemistryModel") >> userModel;

        // construct chemistry model type name by inserting first template
        // argument
        label tempOpen = userModel.find('<');
        label tempClose = userModel.find('>');

        word className = userModel(0, tempOpen);
        thermoTypeName = userModel(tempOpen + 1, tempClose - tempOpen - 1);

        ePsiChemistryModelType =
            className + '<' + typeName + ',' + thermoTypeName + '>';
    }

    if (debug)
    {
        Info<< "Selecting ePsiChemistryModel " << ePsiChemistryModelType << endl;
    }
    else
    {
        Info<< "Selecting ePsiChemistryModel " << userModel << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(ePsiChemistryModelType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        if (debug)
        {
            FatalErrorIn("ePsiChemistryModelBase::New(const mesh&)")
                << "Unknown ePsiChemistryModel type " << ePsiChemistryModelType
                << nl << nl << "Valid ePsiChemistryModel types are:" << nl
                << fvMeshConstructorTablePtr_->sortedToc() << nl << exit(FatalError);
        }
        else
        {
            wordList models = fvMeshConstructorTablePtr_->sortedToc();
            forAll(models, i)
            {
                models[i] = models[i].replace(typeName + ',', "");
            }

            FatalErrorIn("ePsiChemistryModelBase::New(const mesh&)")
                << "Unknown ePsiChemistryModel type " << userModel
                << nl << nl << "Valid ePsiChemistryModel types are:" << nl
                << models << nl << exit(FatalError);
        }
    }

        Info<< "Returning from ePsiChemistryModel::New" << endl;
    return autoPtr<ePsiChemistryModel>
        (cstrIter()(mesh, typeName, thermoTypeName));
}


// ************************************************************************* //
