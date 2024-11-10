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
    Foam::eCombustionThermo

Description
    Foam::eCombustionThermo

SourceFiles
    eCombustionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eCombustionThermo_H
#define eCombustionThermo_H

#include "basicPsiThermo.H"
#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class eCombustionThermo Declaration
\*---------------------------------------------------------------------------*/

class eCombustionThermo
:
    public basicPsiThermo
{

protected:

    // Protected data

        //- Enthalpy field
        volScalarField e_;

        volScalarField es_;
        volScalarField h_;
        volScalarField hs_;



public:

    //- Runtime type information
    TypeName("eCombustionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        eCombustionThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from dictionary and mesh
        eCombustionThermo(const fvMesh&);


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<eCombustionThermo> New(const fvMesh&);

        //- Select and check that package contains 'thermoType'
        static autoPtr<eCombustionThermo> NewType
        (
            const fvMesh&,
            const word& thermoType
        );


    //- Destructor
    virtual ~eCombustionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicMultiComponentMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const = 0;


        // Access to thermodynamic state variables

            //- Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& e()
            {
                return e_;
            }

            //- Internal energy [J/kg]
            virtual const volScalarField& e() const
            {
                return e_;
            }

            virtual volScalarField& hs()
            {  // reimplemented from basicThermo. Required for ODEchemistry
	//Info << "eCombustionThermo: using hs() 1" << endl;
                return hs_;
            }
            virtual volScalarField& es()
            {  // reimplemented from basicThermo. Required only for postprocessing
                return es_;
            }
            virtual volScalarField& h()
            {  // reimplemented from basicThermo. Required only for postprocessing
                return h_;
            }
/*
            virtual volScalarField& hc()
            {  // Required only for postprocessing
                return h_-hs_;
            }
            virtual volScalarField& ec()
            {  // Required only for postprocessing
                return e_-es_;
            }
*/

	    // compute internal energy for a different temperature
            //virtual tmp<volScalarField> e
/*	    
	    tmp<volScalarField> e
            (
                volScalarField& Tgiven
            ) const;
*/
            //- Internal energy for patch [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const label patchi
            ) const;


	    //virtual volScalarField eFromT(const volScalarField& Tgiven) const;	// for MinMax Model: evaluate e at given temperature
	    virtual void eFromT(volScalarField& eToCompute, const volScalarField& Tgiven) const;
	    virtual void hFromT(volScalarField& hToCompute, const volScalarField& Tgiven) const;	    

        //- Update properties
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
