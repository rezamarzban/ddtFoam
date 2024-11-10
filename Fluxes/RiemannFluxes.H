/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    RiemannFlux

Description
    Generic Godunov flux class. Does an extrapolation of cell centered
    variables to face centeres, using a Taylor-Series Expansion. Therefore a
    slope limiter is utilized.

    TODO: Make the slope limiter run-time selectable.

    References:
    Weiss, Jonathan M. and Wayne A. Smith, "Preconditioning Applied to Variable
    and Constant Density Flows." AIAA Journal 33/11 (November 1995): 2050–2057.

Author
    Oliver Borm  All rights reserved.
    Florian Ettner

SourceFiles
    RiemannFlux.H
    RiemannFlux.C

\*---------------------------------------------------------------------------*/

#ifndef RiemannFlux_H
#define RiemannFlux_H

#include "fvCFD.H"
#include "eCombustionThermo.H"
#include "turbulenceModel.H"

//#include "SlopeLimiter/BarthJespersenSlopeMulti/BarthJespersenSlopeMulti.H"
//#include "SlopeLimiter/VenkatakrishnanSlopeMulti/VenkatakrishnanSlopeMulti.H"
//#include "SlopeLimiter/constantSlopeMulti/constantSlopeMulti.H"	// use a constant value instead of a limiter


// #include "NVDTVD.H"
// #include "NVDVTVDV.H"

#include "hllcALEFlux/hllcALEFlux.H"
//#include "roeALEFlux/roeALEFlux.H"
#include "AUSMplusALEFlux/AUSMplusALEFlux.H"
#include "AUSMplusUpALEFlux/AUSMplusUpALEFlux.H"
#include "rusanovALEFlux/rusanovALEFlux.H"

#include "limitedSurfaceInterpolationScheme.H"

//#include "fePtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RiemannFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Flux>
class RiemannFlux
:
    public Flux
{
    // Private data

	    //- Reference to mesh
	    const fvMesh& mesh_;
	    
            //- Reference to the thermophysicalModel
            const eCombustionThermo& thermophysicalModel_;

            //- Reference to the turbulenceModel
            const compressible::turbulenceModel& turbulenceModel_;

	    //- Number of fields
	    const label Npassive_;
	    
            //- References to transported scalars
	    const  volScalarField& rho_;
	    const  volScalarField& rhoE_;
	    const  volScalarField& rhoEu_;
	    
	    const UPtrList<volScalarField>& rhoPassiveScalar_;
	    
            //- Velocity
            const  volVectorField& rhoU_;

            //- static Pressure
	    volScalarField p_; // reqired to compute gradp_
	    
	    // Fluxes
	    surfaceScalarField rhoFlux_;
	    surfaceScalarField rhoEFlux_;
	    surfaceScalarField rhoEuFlux_;	    
	    surfaceScalarField pFlux_;	// rhoHFlux = rhoEFlux_ + pFlux_
	    
	    
            //- Velocity flux
            surfaceVectorField rhoUFlux_;

        // Gradients

            volVectorField gradp_;	    
	    
	    volVectorField gradrho_;
	    volVectorField gradrhoE_;
	    volVectorField gradrhoEu_;
	    
	    
            //- Velocity
            volTensorField gradrhoU_;
	    
    	    //- speed of sound
	    volScalarField a_;	    
            volVectorField grada_;

	    UPtrList<surfaceScalarField> rhoScalarFlux_;	    
	    PtrList<volVectorField> gradrhoScalar_;	

	    // Limiter control
	        // Switch between 1D and MultiDimensional Limiters
	    bool multidimLimiterSwitch;
	    word limiterName; // 1D Limiters
	    word epsilon; // epsilon for the VK limiter
	    scalar Konstant;// constant for entropy fix (Roe) or Ma_Inf (AUSM+up)
	    
	    //volScalarField minimumLimiter;  // make it a member to return a reference for visualization

	  
	    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        RiemannFlux(const RiemannFlux&);

        //- Disallow default bitwise assignment
        void operator=(const RiemannFlux&);



public:

    // Constructors

        //- Construct from components
        RiemannFlux		// rhoU, rho, rhoE and rhoEu need to be put first. The remaining volScalarFields are passive scalars in any order
        (
	  const volVectorField& rhoU,        
	  const volScalarField& rho,
	  const volScalarField& rhoE,
	  const volScalarField& rhoEu,	 
	  const UPtrList<volScalarField>& rhoPassiveScalars,
	  const eCombustionThermo& thermophysicalModel,
	  const compressible::turbulenceModel& turbulenceModel
        );


    // Destructor - default
    
    //  ~RiemannFlux();

    // Member Functions

	    //- Return reference to the mesh
	    const fvMesh& mesh() const
	    {
		return mesh_;
	    }

        // Return fluxes

            const surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            const surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }

            const surfaceScalarField& rhoEuFlux() const
            {
                return rhoEuFlux_;
            }


            const surfaceScalarField& pFlux() const
            {
                return pFlux_;
            }
            const surfaceScalarField& rhoHFlux() const
            {
                return rhoEFlux_+pFlux_;
            }
            

            
            
            const surfaceScalarField& rhoScalarFlux(label i) const
            {
		return rhoScalarFlux_[i];
            }

            //- Return velocity flux
            const surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }
            
            /*
            //- Return velocity flux
            const volScalarField& minimumRiemannLimiter() const
            {
                return minimumLimiter;
            }
            */


        //- Update fluxes based on current state
        void update(Switch secondOrder);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "RiemannFluxes.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
