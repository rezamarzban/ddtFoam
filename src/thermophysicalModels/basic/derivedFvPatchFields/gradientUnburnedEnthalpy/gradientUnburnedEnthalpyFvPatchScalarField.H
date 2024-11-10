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
    Foam::gradientUnburnedEnthalpyFvPatchScalarField

Description
    Gradient boundary condition for enthalpy

SourceFiles
    gradientUnburnedEnthalpyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef gradientUnburnedEnthalpyFvPatchScalarField_H
#define gradientUnburnedEnthalpyFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class gradientUnburnedEnthalpyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class gradientUnburnedEnthalpyFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("gradientUnburnedEnthalpy");


    // Constructors

        //- Construct from patch and internal field
        gradientUnburnedEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gradientUnburnedEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given gradientUnburnedEnthalpyFvPatchScalarField
        // onto a new patch
        gradientUnburnedEnthalpyFvPatchScalarField
        (
            const gradientUnburnedEnthalpyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gradientUnburnedEnthalpyFvPatchScalarField
        (
            const gradientUnburnedEnthalpyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gradientUnburnedEnthalpyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gradientUnburnedEnthalpyFvPatchScalarField
        (
            const gradientUnburnedEnthalpyFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new gradientUnburnedEnthalpyFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
