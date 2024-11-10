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
    Foam::eUnburnedThermo

Description
    copied from basicThermo
    calculates the unburned properties
    unlike the burned thermo class, it requires an up-to-date pressure (to calculate the correct unburned density)
    

SourceFiles
    eUnburnedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eUnburnedThermo_H
#define eUnburnedThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class eUnburnedThermo Declaration
\*---------------------------------------------------------------------------*/

class eUnburnedThermo
:
    public IOdictionary
{

public:

        static const int nCoeffs_=14;
//        typedef scalar coeffArray[nCoeffs_];

        scalar yO2inAir_;

private:

    // Private data
	const scalar MH2_;
    	const scalar MO2_;
    	const scalar MN2_;
	const scalar Runiv_;

        // Temperature limits of applicability of functions
        const scalar Tlow_;
	const scalar Thigh_;
	const scalar Tcommon_;
	const label maxIter_;
        const scalar Ttol_; // Kelvin absolute
        const scalar Tmin_; // bounding
        const scalar Tmax_; // bounding        

	scalar H2coeffs_[nCoeffs_];
	scalar O2coeffs_[nCoeffs_];
	scalar N2coeffs_[nCoeffs_];



    // Private member functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;



protected:

    // Protected data

        // Fields

            //- Temperature [K]
            volScalarField Tu_;

            //- internal energy [J/kg]
            volScalarField eu_;

            //- pressure [N/m^2]
            const volScalarField& p_;  // reference to external field, same pressure in the burned and unburned gas
            
            //- mixture fraction [J/kg]
            const volScalarField& fH_;  // reference to external field
                
            //- gas constant [J/kgK]
            volScalarField Ru_;

            //- compressibility: 1/(Ru*Tu) [kg/J= s^2/m^2]
            volScalarField psiu_;
	    
            //- unburned density [kg/m^3]
            volScalarField rhou_;



    // Protected member functions

        // Enthalpy

            //- Return the unburned internal energy field boundary types by interrogating the
            //  temperature field boundary types
	// not strictly correct, the fuel mixture fraction would also be required
            wordList euBoundaryTypes();

            //- Correct the unburned internal energy field boundaries
            void euBoundaryCorrection(volScalarField& e);


        //- Construct as copy (not implemented)
        eUnburnedThermo(const eUnburnedThermo&);


public:



    //- Runtime type information
    TypeName("eUnburnedThermo");


    // Constructors

        //- Construct from mesh

        eUnburnedThermo(const fvMesh&, const volScalarField &p, const volScalarField &fH, const scalar& yO2inAir);


    //- Destructor
        ~eUnburnedThermo();


    // Member functions

	// Initialize: get hu from Tu
	void initialize();

        //- Update properties: get Tu from hu
        void correct();
	
        //- Update properties: get Tu from hu, while Tu remains bounded by T
	//void correct(const volScalarField& T);

            //- internal energy [J/kg]
            //  Non-const access allowed for transport equations
            volScalarField& eu()
            {
                return eu_;
            }

            //- internal energy [J/kg]
            const volScalarField& hu() const
            {
                return eu_;
            }

        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            const volScalarField& Tu() const
	    {
		return Tu_;
 	    }

            //- Mixture fraction [-]	// function required for hu boundary condition
            const volScalarField& fH() const
            {
                return fH_;
            }
	    
            const volScalarField& psiu() const
	    {
		return psiu_;
 	    }
 	    
            // unburned density
            const volScalarField& rhou() const
	    {
		return rhou_;
	    }
 	    
 	    



            //- Unburned internal energy [J/kg]
	    scalar eu(const scalar& Tui, const scalar& fHi) const;

            //- Unburned internal energy for patch [J/kg]
            tmp<scalarField> eu
            (
                const scalarField& Tu,
                const scalarField& fH,
                const label patchi
            ) const;

            //- Unburned internal energy for cell-set [J/kg]
            tmp<scalarField> eu
            (
                const scalarField& Tu,
                const scalarField& fH,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
	    scalar Cvu(const scalar& Tui, const scalar& fHi) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            tmp<scalarField> Cvu
            (
                const scalarField& Tu,
                const scalarField& fH,
                const label patchi
            ) const;
	    
            //- Heat capacity at constant pressure [J/kg/K]
            tmp<volScalarField> Cvu() const;
	    
	    
	    void euFromT(volScalarField& euToCompute, const volScalarField& Tgiven) const;
	    

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
