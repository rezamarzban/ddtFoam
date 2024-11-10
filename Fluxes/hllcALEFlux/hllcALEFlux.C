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

Application
    hllcALEFlux

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "hllcALEFlux.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //+


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcALEFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    scalar& rhoEuFlux, 
    scalar& pFlux,
    scalarList& rhoScalarFlux,
 
    const scalar& pLeft,
    const scalar& pRight,
 
    const vector& rhoULeft,
    const vector& rhoURight,

    const scalar& rhoLeft,
    const scalar& rhoRight,

    const scalar& aLeft,
    const scalar& aRight,

    const scalar& rhoELeft,
    const scalar& rhoERight,
 
    const scalar& rhoEuLeft,
    const scalar& rhoEuRight,
 
    const scalarList& rhoScalarLeft,
    const scalarList& rhoScalarRight,

    const vector& Sf,
    const scalar& magSf,
    const scalar& K_Roe
) const
{
	//const label N = rhoScalarFlux.size();

        // bounding variables
        const scalar rhoMin = SMALL;

        // Step 1: decode left and right:
        // normal vector
        const vector normalVector = Sf/magSf;

	/*
        // DensityTotalEnergy
        const scalar rhoELeft =  rhoLeft *(eLeft + 0.5*magSqr(ULeft) + kLeft);
        const scalar rhoERight = rhoRight*(eRight + 0.5*magSqr(URight) + kRight);

	scalarList rhoYLeft(species);
	scalarList rhoYRight(species);
	forAll(rhoYFlux,i)
	{
	    rhoYLeft[i]  = rhoLeft*yLeft[i];
	    rhoYRight[i] = rhoRight*yRight[i];    
	}
	    

        // DensityVelocity
        const vector rhoULeft  = rhoLeft *ULeft;
        const vector rhoURight = rhoRight*URight;
*/
	const vector ULeft = rhoULeft/rhoLeft;
	const vector URight = rhoURight/rhoRight;
	
        // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
        const scalar qLeft  = ULeft  & normalVector;
        const scalar qRight = URight & normalVector;

        // Step 2:
        // needs rho_{l,r}, U_{l,r}, H_{l,r}, kappa_{l,r}, Gamma_{l,r}, q_{l,r}
        // compute Roe weights
        const scalar rhoLeftSqrt  = Foam::sqrt(max(rhoLeft ,rhoMin));
        const scalar rhoRightSqrt = Foam::sqrt(max(rhoRight,rhoMin));

        const scalar wLeft = rhoLeftSqrt / max((rhoLeftSqrt + rhoRightSqrt),rhoMin);
        const scalar wRight = 1.0 - wLeft;

        // Roe averaged velocity
        const vector UTilde = wLeft*ULeft + wRight*URight;

        // Roe averaged relative face velocity
        const scalar UNormalTilde = UTilde & normalVector;

        // Roe averaged contravariant velocity
//        const scalar contrUTilde = (UTilde & normalVector);

        // Roe averaged total enthalpy
//        const scalar HTilde = wLeft*HLeft + wRight*HRight;

        // Roe averaged kappa
        // TODO: needs to be verified for non constant kappa values
        //const scalar kappaTilde = wLeft*kappaLeft + wRight*kappaRight;

        // Speed of sound with Roe reconstruction values
        //const scalar aTilde =
        //    Foam::sqrt(max(SMALL,(kappaTilde-1.0)*(HTilde-0.5*sqr(contrUTilde))));

	//const scalar aTilde = wLeft*aLeft + wRight*aRight; // passt das?

        // Step 3: compute signal speeds for face:
        //const scalar SLeft  = min(qLeft -aLeft,  UNormalTilde-aTilde);
        //const scalar SRight = max(qRight+aRight, UNormalTilde+aTilde);
        // See for an alternative computation of wave speed: Toro (2009) Kap. 10.5.1
	
	const scalar eta2 = 0.5*rhoLeftSqrt*rhoRightSqrt*pow(rhoLeftSqrt+rhoRightSqrt,-2);
	const scalar d2 = (rhoLeftSqrt*pow(aLeft,2)+rhoRightSqrt*pow(aRight,2))/(rhoLeftSqrt+rhoRightSqrt) + eta2 * pow(qRight-qLeft,2);
	const scalar d2Sqrt = Foam::sqrt(max(d2 ,rhoMin));
	const scalar SLeft = UNormalTilde-d2Sqrt; // see Einfeldt 1988 (SIAM J Num. Anal.) Eq. 5.2
	const scalar SRight = UNormalTilde+d2Sqrt;
	
	
	/*
	// not recommende (Toro p. 328):
	const scalar SLeft = qLeft - aLeft; 
	const scalar SRight = qRight + aRight;
	*/
	
	
	


        const scalar SStar = (rhoRight*qRight*(SRight-qRight)
            - rhoLeft*qLeft*(SLeft - qLeft) + pLeft - pRight )/
            stabilise((rhoRight*(SRight-qRight)-rhoLeft*(SLeft-qLeft)),SMALL);

        // compute pressure in star region from the right side
        const scalar pStarRight =
            rhoRight * (qRight - SRight) * (qRight - SStar) + pRight;

        // should be equal to the left side
        const scalar pStarLeft  =
            rhoLeft  * (qLeft -  SLeft)  * (qLeft - SStar)  + pLeft;

        // give a warning if this is not the case; the reason is usually inconsistency between primitive (p) and conservative variables (rhoX)
        if ( mag(pStarRight-pStarLeft) > 1e-6 )
        {
                Info << "HLLC warning: mag(pStarRight-pStarLeft) too big: " << mag(pStarRight-pStarLeft) << endl;
		Info << "\t" << rhoLeft << "\t" << qLeft << "\t" << SLeft << "\t" << SStar  << "\t" << pLeft << endl;
		Info << "\t" << rhoRight << "\t" << qRight << "\t" << SRight << "\t" << SStar  << "\t" << pRight << endl;
		
        }

        // in theory, pStarRight == pStarLeft
        const scalar pStar = 0.5*(pStarRight+pStarLeft);

        // Step 4: upwinding - compute states:
        scalar convectionSpeed = 0.0;
        scalar rhoState = 0.0;
        vector rhoUState = vector::zero;
        scalar rhoEState = 0.0;
        scalar rhoEuState = 0.0;
        scalar pState = 0.0;
	scalarList rhoScalarState(rhoScalarFlux.size(),0.0);

        // TODO: Maybe one can use pos/neg implementation, but then one has to
        // evaluate all 4 states at each iteration!
        // label A = pos(SLeft);
        // label B = pos(SStar);
        // label C = pos(SRight);
        // please double check the bool operators again, if one want's to
        // implement this!!!
        // scalar convectionSpeed = A*B*C*qLeft+(1-A)*B*C*SStar
        //     +(1-A)*(1-B)*C*SStar+(1-A)*(1-B)*(1-C)*qRight:

        if ( pos(SLeft) )
        {
            // compute F_l
//	    Info << " F1 ";
            convectionSpeed = qLeft;
            rhoState  = rhoLeft;
            rhoUState = rhoULeft;
            rhoEState = rhoELeft;
            rhoEuState = rhoEuLeft;
            pState = pLeft;
	    forAll(rhoScalarFlux,i)
	    {
	      rhoScalarState[i] = rhoScalarLeft[i];
	    }
	}
        else if ( pos(SStar) )
        {
            scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar),SMALL);

            // compute left star region
//	    Info << " F2 "; //Info << " " << SLeft << " " << SStar << " " << omegaLeft << "\t" ;
	    //Info << " " << rhoScalarLeft[1] << " " << rhoLeft << "\t" ;
            convectionSpeed = SStar;
            rhoState  = omegaLeft *   (SLeft - qLeft) * rhoLeft;
	    //Info << "(" << omegaLeft << " * (" << SLeft << " - " << qLeft << ") * " << rhoLeft << endl;	    
            rhoUState = omegaLeft * ( (SLeft - qLeft) * rhoULeft
                + (pStar - pLeft)*normalVector );
            rhoEState = omegaLeft * ( (SLeft - qLeft) * rhoELeft
                - pLeft*qLeft + pStar*SStar );
            rhoEuState = omegaLeft * ( (SLeft - qLeft) * rhoEuLeft
                - pLeft*qLeft + pStar*SStar );
            pState = pStar;
	    forAll(rhoScalarFlux,i)
	    {
		rhoScalarState[i] = omegaLeft * (SLeft - qLeft) * rhoScalarLeft[i];
	    }
   
	    
        }
        else if ( pos(SRight) )
        {
//	    Info << " F3 " ;
            scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar),SMALL);
	    //Info <<  omegaRight << " = 1.0 / (" << SRight << " - " << SStar << ")"<< endl;
            // compute right star region
            convectionSpeed = SStar;
            rhoState  = omegaRight *   (SRight - qRight) * rhoRight;
	    //Info << "(" << omegaRight << " * (" << SRight << " - " << qRight << ") * " << rhoRight << endl;
            rhoUState = omegaRight * ( (SRight - qRight) * rhoURight
                + (pStar - pRight)*normalVector );
            rhoEState = omegaRight * ( (SRight - qRight) * rhoERight
                - pRight*qRight + pStar*SStar );
            rhoEuState = omegaRight * ( (SRight - qRight) * rhoEuRight
                - pRight*qRight + pStar*SStar );
            pState = pStar;
	    forAll(rhoScalarFlux,i)
	    {
		rhoScalarState[i] = omegaRight * (SRight - qRight) * rhoScalarRight[i];
	    }
        }
        else if ( neg(SRight) )
        {
            // compute F_r
//	    Info << " F4 " ;
            convectionSpeed = qRight;
            rhoState  = rhoRight;
            rhoUState = rhoURight;
            rhoEState = rhoERight;
            rhoEuState = rhoEuRight;
            pState = pRight;
	    forAll(rhoScalarFlux,i)
	    {
	      rhoScalarState[i] = rhoScalarRight[i];
	    }
        }
        else
        {
            Info << "   hllcALEFLUX.C: Error in HLLC Riemann solver" << endl;
        }

        rhoFlux  = (convectionSpeed*rhoState)*magSf;
        rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
        rhoEFlux = (convectionSpeed*rhoEState)*magSf;
        rhoEuFlux = (convectionSpeed*rhoEuState)*magSf;
        pFlux 	 = (convectionSpeed*pState)*magSf;
	//rhoHuFlux= (convectionSpeed*(rhoEuState+pState))*magSf;
	forAll(rhoScalarFlux,i)
	{
	    rhoScalarFlux[i] = (convectionSpeed*rhoScalarState[i])*magSf;
	}
	/*
	if (rhoFlux!=0.0) 
	{
	  Info << " rhoFlux = " << rhoFlux << "\t rhoLeft = " << rhoLeft << "\t rhoRight = " << rhoRight << endl;
	  //if( mag(rhoScalarFlux[1]/rhoFlux-0.001) > 1e-5)
	  {
	      //Info << nl << "   " << convectionSpeed << "\t" << rhoScalarState[1] << "\t" << rhoState << "\t" << magSf << endl;
	      
	  }
	}
	*/
	
	
	//Info << endl;
	//Info << " HLLC: " << convectionSpeed << "\t" << rhoState << "\t" << magSf << "\t = " << rhoFlux << endl;
}

// ************************************************************************* //
