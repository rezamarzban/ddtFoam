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
    RiemannFluxes

Description
    generic Godunov flux class
    Does not include limiters. The user has to limit the gradients in the fvSchemes dictionary, e.g.
	gradSchemes
	{
	  default cellMDLimited Gauss linear 1;
	}
	
Author
    Oliver Borm  All rights reserved.
    Florian Ettner

\*---------------------------------------------------------------------------*/

#include "RiemannFluxes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::RiemannFlux<Flux>::RiemannFlux
(
    const volVectorField& rhoU,
    const volScalarField& rho,
    const volScalarField& rhoE,
    const volScalarField& rhoEu,
    const UPtrList<volScalarField>& rhoPassiveScalars,
    const eCombustionThermo& thermophysicalModel,
    const compressible::turbulenceModel& turbulenceModel
)
:
    mesh_(rho.mesh()),
    thermophysicalModel_(thermophysicalModel),
    turbulenceModel_(turbulenceModel),
    Npassive_(rhoPassiveScalars.size()),
   
    rho_(rho),
    rhoE_(rhoE),
    rhoEu_(rhoEu),
    rhoPassiveScalar_(rhoPassiveScalars),
   
    rhoU_(rhoU),
    
    p_(rho_/thermophysicalModel.psi()),
    
    rhoFlux_
    (
        IOobject
        (
            "rhoFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(rhoU_) & mesh_.Sf()         // only initialisation
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	rhoFlux_*linearInterpolate(rhoE_/rho_)         // only initialisation
    ),
    rhoEuFlux_
    (
        IOobject
        (
            "rhoEuFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	rhoFlux_*linearInterpolate(rhoEu_/rho_)         // only initialisation
    ),
    pFlux_
    (
        IOobject
        (
            "pFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	rhoFlux_*linearInterpolate(p_/rho_)         // only initialisation
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(rhoU_/rho_)         // only initialisation
    ),

    gradp_(fvc::grad(p_,"grad(scalarSlope)")),
    gradrho_(fvc::grad(rho_,"grad(scalarSlope)")),
    gradrhoE_(fvc::grad(rhoE_,"grad(scalarSlope)")), 
    gradrhoEu_(fvc::grad(rhoEu_,"grad(scalarSlope)")), 
    gradrhoU_(fvc::grad(rhoU_,"grad(USlope)")),
    a_(sqrt(thermophysicalModel_.Cp()/thermophysicalModel_.Cv() / thermophysicalModel_.psi())),
    grada_(fvc::grad(a_,"grad(scalarSlope)")) /* ,
    
    minimumLimiter
    (
        IOobject
        (
            "minimumLimiter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    )   
    */
{
      Info << "setting " << endl;
      Info << "Number of passive scalars: " << Npassive_ << endl;

      rhoScalarFlux_.setSize(Npassive_);
      gradrhoScalar_.setSize(Npassive_);
    
  
// these are variables  => & operator required, otherwise a copy is created
      //rhoScalarFlux_.set(0,&rhoCFlux_);	// results in a lost pointer in the destructor
      //rhoScalarFlux_.set(1,&rhoTauFlux_);
      //rhoScalarFlux_.set(2,&rhofHFlux_);

      forAll(rhoPassiveScalar_,i)
      {	    
	    rhoScalarFlux_.set
	    (
	      i, 
	      new surfaceScalarField
	      (
		IOobject
		(
		  "rhoScalarFlux",
		  mesh_.time().timeName(),
		  mesh_,
		  IOobject::NO_READ,
		  IOobject::NO_WRITE
		),
		rhoFlux_*linearInterpolate(rhoPassiveScalar_[i]/rho_)         // only initialisation
	      )
	    );	    
	    
      }	

	    forAll(rhoPassiveScalar_,i)
	    {	
	    gradrhoScalar_.set
	    (
	      i,
	      new volVectorField
	      (
		fvc::grad(rhoPassiveScalar_[i],"grad(scalarSlope)")
	      )
	    );
	    }
// these are variables  => & operator required, otherwise a copy is created
      //gradrhoScalar_.set(0,&gradrhoC_);
      //gradrhoScalar_.set(1,&gradrhoTau_);
      //gradrhoScalar_.set(2,&gradrhofH_);
      
          multidimLimiterSwitch=false;  // only intialisation
	  limiterName="vanAlbadaSlope"; // only intialisation
	  epsilon="5";
	  Konstant=0.05;	// MaInf (AUSM+up) or entropy fix parameter (Roe)
	  
    // read riemann solver coeffs
    if(mesh_.solutionDict().found("Riemann"))   // system/fvSolution
    {
        dictionary riemann = mesh_.solutionDict().subDict("Riemann");
        if (riemann.found("multidimLimiter"))
        {
            multidimLimiterSwitch = Switch(riemann.lookup("multidimLimiter"));
        }
        if (riemann.found("limiterName"))
        {
            limiterName = word(riemann.lookup("limiterName"));
        }
        if (riemann.found("epsilon"))
        {
            epsilon = word(riemann.lookup("epsilon"));
        }
        Konstant = riemann.lookupOrDefault("RoeKonstant",Konstant);
    }


      
      Info << "end of constructor " << endl;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux>
void Foam::RiemannFlux<Flux>::update(Switch secondOrder)
{   
 
    Info << "Calculating Riemann fluxes, second order = " << secondOrder << endl;
    //Info << "rhoE max. = " << max(rhoE_).value() << endl;
    
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    const volVectorField& cellCenter = mesh_.C();
    const surfaceVectorField& faceCenter = mesh_.Cf();
   
    a_ = sqrt(thermophysicalModel_.Cp()/thermophysicalModel_.Cv() / thermophysicalModel_.psi()); // speed of sound
    p_ = rho_/thermophysicalModel_.psi();

    

    //if (secondOrder)
    {
	// compute gradients
	gradrho_  = fvc::grad(rho_,"grad(scalarSlope)");
	
	/*
	volVectorField linearRhoGradient = fvc::grad(rho_,"linear");
	//linearRhoGradient.max(SMALL);
	minimumLimiter = mag(gradrho_) / max(mag(linearRhoGradient),dimensionedScalar("dummy",dimMass/dimVol/dimLength,SMALL));
	minimumLimiter.min(1.0);
	//Info<< " rhoLimiter min...max   = " << min(minimumLimiter).value() << " ... " << max(minimumLimiter).value() << endl;    	
	*/
	gradrhoE_ = fvc::grad(rhoE_,"grad(scalarSlope)");
	gradrhoEu_ = fvc::grad(rhoEu_,"grad(scalarSlope)");
	
	
	grada_ = fvc::grad(a_,"grad(scalarSlope)");
	gradp_ = fvc::grad(p_,"grad(scalarSlope)");
	
	gradrhoU_ = fvc::grad(rhoU_,"grad(USlope)");  
	
	forAll(rhoPassiveScalar_,i) 
	{
	    gradrhoScalar_[i] = fvc::grad(rhoPassiveScalar_[i],"grad(scalarSlope)");	
	}
      
    } 


    scalarList rhoScalarOwns(Npassive_);
    scalarList rhoScalarNeis(Npassive_);
    scalarList rhoScalarFluxList(Npassive_);
   

  
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector deltaRLeft  = faceCenter[faceI] - cellCenter[own];
        vector deltaRRight = faceCenter[faceI] - cellCenter[nei];
	
	forAll(rhoScalarOwns,i)
	{
	    const volScalarField& rhoScalari = rhoPassiveScalar_[i];
	    const volVectorField& gradrhoScalari = gradrhoScalar_[i];

	    rhoScalarOwns[i] = rhoScalari[own] + secondOrder*(deltaRLeft  & gradrhoScalari[own]); 
	    rhoScalarNeis[i] = rhoScalari[nei] + secondOrder*(deltaRRight & gradrhoScalari[nei]); 
	}

        // calculate fluxes with reconstructed conservative variables at faces 
        Flux::evaluateFlux // compute new values of rhoUFlux_ and rhoScalarFluxList
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            rhoEuFlux_[faceI],
            pFlux_[faceI],
	    rhoScalarFluxList,
	 
	    p_[own] + secondOrder*(deltaRLeft  & gradp_[own]),
            p_[nei] + secondOrder*(deltaRRight & gradp_[nei]),   
	 
	    rhoU_[own] + secondOrder*(deltaRLeft & gradrhoU_[own]), // bugfix Saegeler delta&U, not U&delta
            rhoU_[nei] + secondOrder*(deltaRRight& gradrhoU_[nei]),

	    rho_[own] + secondOrder*(deltaRLeft  & gradrho_[own]), 
            rho_[nei] + secondOrder*(deltaRRight & gradrho_[nei]),
      
	    a_[own] + secondOrder*(deltaRLeft  & grada_[own]),       
            a_[nei] + secondOrder*(deltaRRight & grada_[nei]),      

	    rhoE_[own] + secondOrder*(deltaRLeft  & gradrhoE_[own]),	
	    rhoE_[nei] + secondOrder*(deltaRRight & gradrhoE_[nei]),
	
	    rhoEu_[own] + secondOrder*(deltaRLeft  & gradrhoEu_[own]),	
	    rhoEu_[nei] + secondOrder*(deltaRRight & gradrhoEu_[nei]),
	 
	    rhoScalarOwns,
	    rhoScalarNeis,

            Sf[faceI],      // face vector
            magSf[faceI],   // face area
            Konstant	    // Roe Konstant (from dictionary, default 0.05)
        );
	
	forAll(rhoScalarOwns,i)	// write results back to rhoScalarFlux_ fields
	{
	    surfaceScalarField& rhoScalarFluxi = rhoScalarFlux_[i];
	    rhoScalarFluxi[faceI] = rhoScalarFluxList[i];
	}
    }
    
   
    // Update boundary field and values
    forAll(rho_.boundaryField(), patchi)
    {      
        fvsPatchScalarField& patchrhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& patchrhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& patchrhoEFlux = rhoEFlux_.boundaryField()[patchi];
        fvsPatchScalarField& patchrhoEuFlux = rhoEuFlux_.boundaryField()[patchi];
        fvsPatchScalarField& patchpFlux    = pFlux_.boundaryField()[patchi];	

        const fvPatchScalarField& patchp = p_.boundaryField()[patchi];
        const fvPatchVectorField& patchrhoU = rhoU_.boundaryField()[patchi];
        const fvPatchScalarField& patchrho = rho_.boundaryField()[patchi];
        const fvPatchScalarField& patcha = a_.boundaryField()[patchi];
	const fvPatchScalarField& patchrhoE = rhoE_.boundaryField()[patchi];	
	const fvPatchScalarField& patchrhoEu = rhoEu_.boundaryField()[patchi];	
	

        const fvPatchVectorField& patchGradp = gradp_.boundaryField()[patchi];
        const fvPatchTensorField& patchGradrhoU = gradrhoU_.boundaryField()[patchi];
        const fvPatchVectorField& patchGradrho = gradrho_.boundaryField()[patchi];
	const fvPatchVectorField& patchGrada = grada_.boundaryField()[patchi];
	const fvPatchVectorField& patchGradrhoE = gradrhoE_.boundaryField()[patchi];
	const fvPatchVectorField& patchGradrhoEu = gradrhoEu_.boundaryField()[patchi];
	

        const fvsPatchVectorField& patchSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& patchMagSf = magSf.boundaryField()[patchi];


        const fvPatchVectorField& patchCellCenter = cellCenter.boundaryField()[patchi];
	
        if (patchrho.coupled())
        {	//Info << "Patch " << patchi << " is coupled." << endl;
            const scalarField patchpLeft  = patchp.patchInternalField();
            const scalarField patchpRight = patchp.patchNeighbourField();

            const vectorField patchrhoULeft  = patchrhoU.patchInternalField();
            const vectorField patchrhoURight = patchrhoU.patchNeighbourField();

            const scalarField patchrhoLeft  = patchrho.patchInternalField();
            const scalarField patchrhoRight = patchrho.patchNeighbourField();

            const scalarField patchaLeft  = patcha.patchInternalField();
            const scalarField patchaRight = patcha.patchNeighbourField();

            const scalarField patchrhoELeft  = patchrhoE.patchInternalField();
            const scalarField patchrhoERight = patchrhoE.patchNeighbourField();
	    
            const scalarField patchrhoEuLeft  = patchrhoEu.patchInternalField();
            const scalarField patchrhoEuRight = patchrhoEu.patchNeighbourField();

            // cell gradients
            const vectorField patchGradpLeft  = patchGradp.patchInternalField();
            const vectorField patchGradpRight = patchGradp.patchNeighbourField();

            const tensorField patchGradrhoULeft  = patchGradrhoU.patchInternalField();
            const tensorField patchGradrhoURight = patchGradrhoU.patchNeighbourField();

            const vectorField patchGradrhoLeft  = patchGradrho.patchInternalField();
            const vectorField patchGradrhoRight = patchGradrho.patchNeighbourField();

	    const vectorField patchGradaLeft  = patchGrada.patchInternalField();
            const vectorField patchGradaRight = patchGrada.patchNeighbourField();
	    
            const vectorField patchGradrhoELeft  = patchGradrhoE.patchInternalField();
            const vectorField patchGradrhoERight = patchGradrhoE.patchNeighbourField();
	    
            const vectorField patchGradrhoEuLeft  = patchGradrhoEu.patchInternalField();
            const vectorField patchGradrhoEuRight = patchGradrhoEu.patchNeighbourField();
	    
	    
            // cell and face centers
            const vectorField faceCenter =  patchrho.patch().Cf();
            const vectorField patchCellCenterLeft  =  patchCellCenter.patchInternalField();
            const vectorField patchCellCenterRight =  patchCellCenter.patchNeighbourField();

            forAll(patchrho, facei)
            {
		//Info << " face " << facei << endl;	      
                vector deltaRLeft  = faceCenter[facei] - patchCellCenterLeft[facei];
                vector deltaRRight = faceCenter[facei] - patchCellCenterRight[facei];
		
		scalarList patchrhoScalarOwns(Npassive_);
		scalarList patchrhoScalarNeis(Npassive_);
		scalarList patchrhoScalarFluxList(Npassive_);
		
		forAll(patchrhoScalarOwns,i)
		{
		      const fvPatchScalarField& patchrhoScalari = rhoPassiveScalar_[i].boundaryField()[patchi];
		      const scalarField patchrhoScalariLeft = patchrhoScalari.patchInternalField();	// references don't work here (first value of scalarField& patchrhoScalariLeft is always 0) - why ??? patchInternalField() returns a temporary field (?)
		      const scalarField patchrhoScalariRight= patchrhoScalari.patchNeighbourField();
	      
		      const fvPatchVectorField& patchGradrhoScalari = gradrhoScalar_[i].boundaryField()[patchi];
		      const vectorField patchGradrhoScalariLeft  = patchGradrhoScalari.patchInternalField();
		      const vectorField patchGradrhoScalariRight = patchGradrhoScalari.patchNeighbourField();

		      patchrhoScalarOwns[i]=patchrhoScalariLeft[facei]  + secondOrder * (patchGradrhoScalariLeft[facei] & deltaRLeft);
		      patchrhoScalarNeis[i]=patchrhoScalariRight[facei] + secondOrder * (patchGradrhoScalariRight[facei] & deltaRRight);

		      //fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		      //patchrhoScalarFluxList[i] = rhoScalariFluxOnBoundary[facei];	      
		}
	
                // Calculate fluxes at coupled boundary faces
                Flux::evaluateFlux
                (
                    patchrhoFlux[facei],
                    patchrhoUFlux[facei],
                    patchrhoEFlux[facei],
                    patchrhoEuFlux[facei],
                    patchpFlux[facei],
	 	    patchrhoScalarFluxList,

		    patchpLeft[facei]  + secondOrder*(patchGradpLeft[facei] & deltaRLeft),         
		    patchpRight[facei] + secondOrder*(patchGradpRight[facei] & deltaRRight),      
		 
                    patchrhoULeft[facei]  + secondOrder*(deltaRLeft & patchGradrhoULeft[facei]), 
                    patchrhoURight[facei] + secondOrder*(deltaRRight & patchGradrhoURight[facei]),
		 
                    patchrhoLeft[facei]  + secondOrder*(patchGradrhoLeft[facei]  & deltaRLeft),        
                    patchrhoRight[facei] + secondOrder*(patchGradrhoRight[facei] & deltaRRight),       
		 
                    patchaLeft[facei]  + secondOrder*(patchGradaLeft [facei] & deltaRLeft),       
                    patchaRight[facei] + secondOrder*(patchGradaRight[facei] & deltaRRight),      
		 
                    patchrhoELeft[facei]  + secondOrder*(patchGradrhoELeft [facei] & deltaRLeft), 
                    patchrhoERight[facei] + secondOrder*(patchGradrhoERight[facei] & deltaRRight),		 

                    patchrhoEuLeft[facei]  + secondOrder*(patchGradrhoEuLeft [facei] & deltaRLeft), 
                    patchrhoEuRight[facei] + secondOrder*(patchGradrhoEuRight[facei] & deltaRRight),		 
		 
		    patchrhoScalarOwns,
		    patchrhoScalarNeis,
		 
                    patchSf[facei],       // face vector
                    patchMagSf[facei],    // face area
                    Konstant
                );
		
		forAll(rhoScalarOwns,i)  // writing the resulting flux from the temporary to the permanent field
		{
		    fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		    rhoScalariFluxOnBoundary[facei] = patchrhoScalarFluxList[i];
		}
	
            }
        }
        else  // patch is not coupled
        {
		scalarList patchrhoScalarOwns(Npassive_);
		scalarList patchrhoScalarNeis(Npassive_);
		scalarList patchrhoScalarFluxList(Npassive_);

            forAll(patchrho, facei)
            {

		forAll(patchrhoScalarOwns,i)
		{
		      const fvPatchScalarField& rhoScalari = rhoPassiveScalar_[i].boundaryField()[patchi];
		      patchrhoScalarOwns[i] = rhoScalari[facei];
		}
		// calculate fluxes at not-coupled boundaries
                Flux::evaluateFlux
                (
                    patchrhoFlux[facei],
                    patchrhoUFlux[facei],
                    patchrhoEFlux[facei],
                    patchrhoEuFlux[facei],
                    patchpFlux[facei],
		    patchrhoScalarFluxList,

		    patchp[facei],        
                    patchp[facei],       
		 
                    patchrhoU[facei],     
                    patchrhoU[facei],     

		    patchrho[facei],      
                    patchrho[facei],       

                    patcha[facei],      
                    patcha[facei],      

                    patchrhoE[facei],
                    patchrhoE[facei],
		 
                    patchrhoEu[facei],
                    patchrhoEu[facei],
		 
		    patchrhoScalarOwns,		// scalar list containing the values for all species
		    patchrhoScalarOwns,		// scalar list containing the values for all species

		    patchSf[facei],       
                    patchMagSf[facei],    
                    Konstant
                );
		/*
		Info << facei << "\t rhoFlux = " << patchrhoFlux << endl;
		Info << facei << "\t rhoUFlux = " << patchrhoUFlux << endl;
		Info << facei << "\t rhoEFlux = " << patchrhoEFlux << endl;
		*/
		forAll(rhoScalarOwns,i)  // write results from rhoScalarFluxList back to the correct position in rhoScalarFlux_
		{
	  	  fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		  rhoScalariFluxOnBoundary[facei] = patchrhoScalarFluxList[i];
		}
		
		//Info << "rhofHFlux on boundary " << patchi << ": " << rhoScalarFlux_[2].boundaryField()[patchi] << endl;
		
	    }
        }
    }
	

 
}

// ************************************************************************* //
