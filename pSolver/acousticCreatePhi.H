/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Global
    acousticCreatePhi

Description
    Creates and initialises the face-flux field acousticPhi computed from U+c

\*---------------------------------------------------------------------------*/

#ifndef acousticCreatePhi_H
#define acousticCreatePhi_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading/calculating face flux field acousticPhi\n" << endl;

volVectorField uPlusC = U;
forAll(uPlusC,i)
{
   for(label j=0; j<3; j++)
   {
     uPlusC[i].component(j) = max(U[i].component(j)+cSound[i],U[i].component(j)-cSound[i]);
   }
}

surfaceScalarField acousticPhi
(
    IOobject
    (
        "acousticPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho*uPlusC) & mesh.Sf()
);


#endif

// ************************************************************************* //
