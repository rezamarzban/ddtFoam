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
    Foam::unburnedThermo

Description
    copied from basicThermo
    calculates the unburned properties
    unlike the burned thermo class, it requires an up-to-date pressure (to calculate the correct unburned density)
    

SourceFiles
    unburnedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef unburnedThermo_H
#define unburnedThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class unburnedThermo Declaration
\*---------------------------------------------------------------------------*/

class unburnedThermo
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

            //- enthalpy [J/kg]
            volScalarField hu_;

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

            //- Return the unburned enthalpy field boundary types by interrogating the
            //  temperature field boundary types
	// not strictly correct, the fuel mixture fraction would also be required
            wordList huBoundaryTypes();

            //- Correct the unburned enthalpy field boundaries
            void huBoundaryCorrection(volScalarField& h);


        //- Construct as copy (not implemented)
        unburnedThermo(const unburnedThermo&);


public:



    //- Runtime type information
    TypeName("unburnedThermo");


    // Constructors

        //- Construct from mesh
//        unburnedThermo(const fvMesh&, const volScalarField&);
        unburnedThermo(const fvMesh&, const volScalarField &p, const volScalarField &fH, const scalar& yO2inAir);


    //- Destructor
        ~unburnedThermo();


    // Member functions

	// Initialize: get hu from Tu
	void initialize();

        //- Update properties: get Tu from hu
        void correct();
	
        //- Update properties: get Tu from hu, while Tu remains bounded by T
	void correct(const volScalarField& T);

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            volScalarField& hu()
            {
                return hu_;
            }

            //- Enthalpy [J/kg]
            const volScalarField& hu() const
            {
                return hu_;
            }

        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            const volScalarField& Tu() const
	    {
		return Tu_;
 	    }

            //- Mixture fraction [-]
            //  Non-const access allowed for transport equations
	    /*
            volScalarField& fH()
            {
                return fH_;
            }
            */

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
 	    
 	    



            //- Unburned enthalpy [J/kg]
	    scalar hu(const scalar& Tui, const scalar& fHi) const;

            //- Unburned enthalpy for patch [J/kg]
            tmp<scalarField> hu
            (
                const scalarField& Tu,
                const scalarField& fH,
                const label patchi
            ) const;

            //- Unburned enthalpy for cell-set [J/kg]
            tmp<scalarField> hu
            (
                const scalarField& Tu,
                const scalarField& fH,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
	    scalar Cpu(const scalar& Tui, const scalar& fHi) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            tmp<scalarField> Cpu
            (
                const scalarField& Tu,
                const scalarField& fH,
                const label patchi
            ) const;
	    
            //- Heat capacity at constant pressure [J/kg/K]
            tmp<volScalarField> Cpu() const;
	    
	    
	    void huFromT(volScalarField& huToCompute, const volScalarField& Tgiven) const;
	    




/*

        // Access to thermodynamic state variables


 
            //- Total enthalpy for cell-set [J/kg]
            tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            const volScalarField& T() const;


        // Access to transport state variables


        //- Read thermophysicalProperties dictionary
        virtual bool read();

*/


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
