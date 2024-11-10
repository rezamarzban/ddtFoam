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
    createPhiC

Description
    Creates and initialises the face-flux field phiC.
    This is required for the correct determination of grad(c) in combustion models
	such as burning velocity models.
    In these cases, the gradient has to be inside the unburned fluid,
	not spread between burned and unburned,
	in order to gain a correct burning velocity.
    Specify in the system/fvSchemes dictionary:
	grad(c)		Gauss upwind phiC

\*---------------------------------------------------------------------------*/

#ifndef createPhiC_H
#define createPhiC_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info << "Creating face flux field phiC" << endl;

volScalarField ctemp = 1.0*c; // required to define a separate grad scheme for ctemp

volVectorField cgradtemp = fvc::grad(ctemp);	// use linear grad scheme, just to get the direction of the gradient

surfaceScalarField phiC	// required for upwind scheme: 		grad(c) Gauss upwind phiC;
(
    IOobject
    (
        "phiC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(-cgradtemp) & mesh.Sf()
);
phiC.write();


#endif

// ************************************************************************* //
