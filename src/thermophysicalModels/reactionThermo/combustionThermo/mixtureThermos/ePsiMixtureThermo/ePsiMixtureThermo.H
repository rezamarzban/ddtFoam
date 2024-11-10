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
    Foam::ePsiMixtureThermo

Description
    Foam::ePsiMixtureThermo

SourceFiles
    ePsiMixtureThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ePsiMixtureThermo_H
#define ePsiMixtureThermo_H

#include "eCombustionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ePsiMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class ePsiMixtureThermo
:
    public eCombustionThermo,
    public MixtureType
{
    // Private member functions

        void calculate();

        //- Construct as copy (not implemented)
        ePsiMixtureThermo(const ePsiMixtureThermo<MixtureType>&);


public:

    //- Runtime type information
    TypeName("ePsiMixtureThermo");


    // Constructors

        //- Construct from mesh
        ePsiMixtureThermo(const fvMesh&);


    //- Destructor
    virtual ~ePsiMixtureThermo();


    // Member functions

        //- Return the compostion of the multi-component mixture
        virtual basicMultiComponentMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();

        //- Chemical internal energy [J/kg]
        //virtual tmp<volScalarField> ec() const;

        //- Chemical enthalpy [J/kg]
        virtual tmp<volScalarField> hc() const;

        //- Sensible enthalpy [J/kg]
//        virtual tmp<volScalarField> hs() const;
//virtual const Foam::volScalarField& hs() const;
            //virtual tmp<volScalarField> hs() const;
//            virtual volScalarField& hs();
//Foam::volScalarField& hs();


        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const labelList& cells
            ) const;
/*	    
	    // Calculate energy for the actual composition at a different temperature
            virtual tmp<volScalarField> e
            (
                volScalarField& Tgiven
            ) const;
*/	    
	    
            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;
	    
	    //virtual volScalarField eFromT(const volScalarField& Tgiven) const;	// reimplemented from eCombustionThermo
	    virtual void eFromT(volScalarField& eToCompute, const volScalarField& Tgiven) const;	// reimplemented from eCombustionThermo
	    virtual void hFromT(volScalarField& hToCompute, const volScalarField& Tgiven) const;	// reimplemented from eCombustionThermo	    


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ePsiMixtureThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
