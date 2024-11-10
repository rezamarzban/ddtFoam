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

Class
    Foam::ePsiChemistryModel

Description
    Chemistry model for compressibility-based thermodynamics

SourceFiles
    ePsiChemistryModelI.H
    ePsiChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ePsiChemistryModel_H
#define ePsiChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "eCombustionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class ePsiChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class ePsiChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        ePsiChemistryModel(const ePsiChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const ePsiChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<eCombustionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("ePsiChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        ePsiChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        ),
        (mesh, compTypeName, thermoTypeName)
    );


    // Constructors

        //- Construct from mesh
        ePsiChemistryModel(const fvMesh& mesh, const word& thermoTypeName);


    //- Selector
    static autoPtr<ePsiChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~ePsiChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline eCombustionThermo& thermo();

        //- Return const access to the thermo package
        inline const eCombustionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ePsiChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
