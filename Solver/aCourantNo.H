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

Global
    aCourantNo

Description
    Calculates the mean and maximum wave speed based Courant Numbers (u+cSound,u-cSound),
	taken from rhoCentralFoam.

\*---------------------------------------------------------------------------*/

	cSound = sqrt(thermo.Cp()/thermo.Cv()/thermo.psi()); 

        dimensionedScalar v_zero("v_zero",dimVolume/dimTime, 0.0); 
	
        surfaceScalarField rho_pos =
            fvc::interpolate(rho, pos, "reconstruct(rho)");
        surfaceScalarField rho_neg =
            fvc::interpolate(rho, neg, "reconstruct(rho)");

        surfaceVectorField rhoU_pos =
            fvc::interpolate(rhoU, pos, "reconstruct(U)");
        surfaceVectorField rhoU_neg =
            fvc::interpolate(rhoU, neg, "reconstruct(U)");

        surfaceVectorField U_pos = rhoU_pos/rho_pos;
        surfaceVectorField U_neg = rhoU_neg/rho_neg;

        surfaceScalarField phiv_pos = U_pos & mesh.Sf();
        surfaceScalarField phiv_neg = U_neg & mesh.Sf();

        surfaceScalarField cSf_pos = fvc::interpolate(cSound, pos, "reconstruct(cSound)")*mesh.magSf();
        surfaceScalarField cSf_neg = fvc::interpolate(cSound, neg, "reconstruct(cSound)")*mesh.magSf();

        surfaceScalarField ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
        surfaceScalarField am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

scalar acousticCoNum = 0.0;
scalar acousticmeanCoNum = 0.0;

if (mesh.nInternalFaces())
{
    surfaceScalarField amaxSfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*amaxSf;

    acousticCoNum = max(amaxSfbyDelta/mesh.magSf()).value()*runTime.deltaT().value();

    acousticmeanCoNum =
        (sum(amaxSfbyDelta)/sum(mesh.magSf())).value()
       *runTime.deltaT().value();
}

Info<< "Acoustic Courant Number mean: " << acousticmeanCoNum 
	<< " max: " << acousticCoNum << endl;

scalar CoNum = acousticCoNum;

