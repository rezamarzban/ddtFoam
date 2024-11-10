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

Class
    Foam::xCombustionThermo

Description
    Foam::xCombustionThermo

SourceFiles
    xCombustionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef xCombustionThermo_H
#define xCombustionThermo_H

#include "basicPsiThermo.H"
#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class xCombustionThermo Declaration
\*---------------------------------------------------------------------------*/

class xCombustionThermo
:
    public basicPsiThermo
{

protected:

    // Protected data

        //- Enthalpy field
        volScalarField h_;


public:

    //- Runtime type information
    TypeName("xCombustionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        xCombustionThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from dictionary and mesh
        xCombustionThermo(const fvMesh&);


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<xCombustionThermo> New(const fvMesh&);

        //- Select and check that package contains 'thermoType'
        static autoPtr<xCombustionThermo> NewType
        (
            const fvMesh&,
            const word& thermoType
        );


    //- Destructor
    virtual ~xCombustionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicMultiComponentMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const = 0;


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }


        //- Update properties
        virtual void correct() = 0;

	void hello();

	//- initialize hunburned
//	virtual void initializeH();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
