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

Class
    Foam::hPsiChemistryModel

Description
    FE: Chemistry model for compressibility-based thermodynamics, based on h, not hs

SourceFiles
    hPsiChemistryModelI.H
    hPsiChemistryModel.C
    newhPsiChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef hPsiChemistryModel_H
#define hPsiChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hCombustionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class hPsiChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class hPsiChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        hPsiChemistryModel(const hPsiChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const hPsiChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<hCombustionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("hPsiChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        hPsiChemistryModel,
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
        hPsiChemistryModel(const fvMesh& mesh, const word& thermoTypeName);


    //- Selector
    static autoPtr<hPsiChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~hPsiChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline hCombustionThermo& thermo();

        //- Return const access to the thermo package
        inline const hCombustionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hPsiChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
