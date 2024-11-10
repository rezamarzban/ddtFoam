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

\*---------------------------------------------------------------------------*/

#include "ePsiMixtureThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::ePsiMixtureThermo<MixtureType>::calculate()
{
//	Info << "Updating ePsiMixture:" << nl 
//	<< " e   min...max   = " << min(e_).value()  << " ... " << max(e_).value() << endl;

// typedef sutherlandTransport<specieThermo<janafThermo<perfectGas> > > gasThermoPhysics

    const scalarField& eCells = e_.internalField();
    const scalarField& pCells = p_.internalField();
    scalarField& TCells = T_.internalField();
    scalarField& psiCells = psi_.internalField();
    scalarField& muCells = mu_.internalField();
    scalarField& alphaCells = alpha_.internalField();

    scalarField& hCells  = h_.internalField();
//    scalarField& esCells = es_.internalField();
    scalarField& hsCells = hs_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture =
            this->cellMixture(celli);		// multiComponentMixture.C: sum up mixture properties over all species mass fractions
	    
        TCells[celli] = mixture.TE(eCells[celli], TCells[celli]);  // iteratively from JANAF polynom
      //Info << " calling psi " << endl;
        psiCells[celli] = mixture.psi(pCells[celli], TCells[celli]); // perfectGasI.H: returns 1.0/(R()*T), p is not used
								     // R() references specie::R()
        muCells[celli] = mixture.mu(TCells[celli]);
        alphaCells[celli] = mixture.alpha(TCells[celli]);

//	esCells[celli] = eCells[celli] - mixture.e(specie::Tstd);
//	esCells[celli] = mixture.Es(TCells[celli]);
//	hCells[celli] = eCells[celli] + 1.0/psiCells[celli];
	hCells[celli] = mixture.H(TCells[celli]);
	hsCells[celli] = mixture.Hs(TCells[celli]);

      /*
      if(mixture.W()<10.0)  // if the molecular weight becomes this low, something seems to go wrong
      {
	  Info << " ePsiMixtureThermo.C: W = " << mixture.W() << " in cell " << celli << endl;
	  //Info << "     " << mixture.speciesData() << endl;
	  //Info << "   " << mixture.Y() << endl;	  
      }
      */
	
    }


    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = p_.boundaryField()[patchi];
        fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = psi_.boundaryField()[patchi];

        fvPatchScalarField& pe = e_.boundaryField()[patchi];

        fvPatchScalarField& ph = h_.boundaryField()[patchi];
        //fvPatchScalarField& pes = es_.boundaryField()[patchi];
        fvPatchScalarField& phs = hs_.boundaryField()[patchi];


        fvPatchScalarField& pmu_ = mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha_ = alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                pe[facei] = mixture.E(pT[facei]);	

                ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture.mu(pT[facei]);
                palpha_[facei] = mixture.alpha(pT[facei]);

		//ph[facei] = pe[facei] + 1.0/ppsi[facei];
		ph[facei]  = mixture.H(pT[facei]);
		//pes[facei] = mixture.Es(pT[facei]);
		phs[facei] = mixture.Hs(pT[facei]);

            }
/*
Info    << "patch " << patchi << ": " << nl
	<< " e  min...max   = " << min(pe) << " ... " << max(pe) << nl
	<< " T  min...max   = " << min(pT) << " ... " << max(pT) << nl
	<< " h  min...max   = " << min(ph) << " ... " << max(ph) << nl
	;
*/
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture.TE(pe[facei], pT[facei]);

                ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture.mu(pT[facei]);
                palpha_[facei] = mixture.alpha(pT[facei]);

		ph[facei]  = mixture.H(pT[facei]);
		//pes[facei] = mixture.Es(pT[facei]);
		phs[facei] = mixture.Hs(pT[facei]);

            }
/*
Info    << "patch " << patchi << ": " << nl
	<< " e  min...max   = " << min(pe) << " ... " << max(pe) << nl
	<< " T  min...max   = " << min(pT) << " ... " << max(pT) << nl
	<< " h  min...max   = " << min(ph) << " ... " << max(ph) << nl
	;
*/
        }
    }

/*
Info    << "internal field: " << nl
	<< " T   min...max   = " << min(T_).value()  << " ... " << max(T_).value() << nl
	<< " h   min...max   = " << min(h_).value()  << " ... " << max(h_).value() << nl
	<< " es  min...max   = " << min(es_).value() << " ... " << max(es_).value() << nl
	<< " hs  min...max   = " << min(hs_).value() << " ... " << max(hs_).value() << nl
	<< "*" << endl
	;
*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::ePsiMixtureThermo<MixtureType>::ePsiMixtureThermo(const fvMesh& mesh)
:
    eCombustionThermo(mesh),
    MixtureType(*this, mesh)
{
    scalarField& eCells = e_.internalField();
    const scalarField& TCells = T_.internalField();

    forAll(eCells, celli)
    {
        eCells[celli] = this->cellMixture(celli).E(TCells[celli]);	// defined in specieThermoI.H
	//Info << " ePsiMixtureThermo: e(" << TCells[celli] << ") = " << eCells[celli] << endl;
    }

    forAll(e_.boundaryField(), patchi)
    {
        e_.boundaryField()[patchi] == e(T_.boundaryField()[patchi], patchi);
    }

    eBoundaryCorrection(e_);

    calculate();

    // Switch on saving old time
    psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::ePsiMixtureThermo<MixtureType>::~ePsiMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::ePsiMixtureThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering ePsiMixtureThermo<MixtureType>::correct()" << endl;
    }

    // force the saving of the old-time values
    psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting ePsiMixtureThermo<MixtureType>::correct()" << endl;
    }
}

/*
template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::ePsiMixtureThermo<MixtureType>::ec() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tec
    (
        new volScalarField
        (
            IOobject
            (
                "ec",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            e_.dimensions()
        )
    );

    volScalarField& ecf = tec();
    scalarField& ecCells = ecf.internalField();

    forAll(ecCells, celli)
    {
        ecCells[celli] = this->cellMixture(celli).Ec();
    }

    forAll(ecf.boundaryField(), patchi)
    {
        scalarField& ecp = ecf.boundaryField()[patchi];

        forAll(ecp, facei)
        {
            ecp[facei] = this->patchFaceMixture(patchi, facei).Ec();
        }
    }

    return tec;
}
*/

template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::ePsiMixtureThermo<MixtureType>::hc() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            e_.dimensions()
        )
    );

    volScalarField& hcf = thc();
    scalarField& hcCells = hcf.internalField();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    forAll(hcf.boundaryField(), patchi)
    {
        scalarField& hcp = hcf.boundaryField()[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }
//Info << "ePsiMixtureThermo: Using hc " << endl;
    return thc;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::ePsiMixtureThermo<MixtureType>::e
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> te(new scalarField(T.size()));
    scalarField& e = te();

    forAll(T, celli)
    {
        e[celli] = this->cellMixture(cells[celli]).E(T[celli]);
    }

    return te;
}

/*
// Calculate energy for the actual composition at a different temperature
template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::ePsiMixtureThermo<MixtureType>::e
(
    volScalarField& Tgiven
) const
{
  Info << "using ePsiMixtureThermo<MixtureType>::e(volScalarField& Tgiven) const" << endl;
    //tmp<volscalarField> te(new volScalarField(Tgiven.size()));
    tmp<volScalarField> te(e_);
    volScalarField& e = te();
  if(Tgiven.size() == e_.size())
  {
    forAll(Tgiven, celli)
    {
        e[celli] = this->cellMixture(celli).E(Tgiven[celli]);
    }
  }
  else
  {
               FatalErrorIn
            (
                "ePsiMixtureThermo<MixtureType>::e "
                "( const scalarField& Tgiven ) const "
            )   << "size of energy field and temperature field provided are not the same"
                << abort(FatalError);
  }
      return te;
}
*/



template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::ePsiMixtureThermo<MixtureType>::e
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> te(new scalarField(T.size()));
    scalarField& e = te();

    forAll(T, facei)
    {
        e[facei] = this->patchFaceMixture(patchi, facei).E(T[facei]);
    }

    return te;
}



template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::ePsiMixtureThermo<MixtureType>::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> ths(new scalarField(T.size()));
    scalarField& hs = ths();

    forAll(T, celli)
    {
        hs[celli] = this->cellMixture(cells[celli]).Hs(T[celli]);
    }

    return ths;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::ePsiMixtureThermo<MixtureType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> ths(new scalarField(T.size()));
    scalarField& hs = ths();

    forAll(T, facei)
    {
        hs[facei] = this->patchFaceMixture(patchi, facei).Hs(T[facei]);
    }

    return ths;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::ePsiMixtureThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));

    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
    }

    return tCp;
}



template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::ePsiMixtureThermo<MixtureType>::Cp() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();

    scalarField& cpCells = cp.internalField();
    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        cpCells[celli] = this->cellMixture(celli).Cp(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
    }

    return tCp;
}

template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::ePsiMixtureThermo<MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));

    scalarField& cv = tCv();

    forAll(T, facei)
    {
        cv[facei] = this->patchFaceMixture(patchi, facei).Cv(T[facei]);
    }

    return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::ePsiMixtureThermo<MixtureType>::Cv() const
{
//  Info << "Calling ePsiMixtureThermo<MixtureType>::Cv() const" << endl;
    const fvMesh& mesh = T_.mesh();
//    Info << "T_ = " << T_ << endl;		// kaputt!
//Info << " OK 01 " << endl;
    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );
    
//Info << " OK 05 " << endl;
    volScalarField& cv = tCv();

    scalarField& cvCells = cv.internalField();

    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        cvCells[celli] = this->cellMixture(celli).Cv(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] = Cv(T_.boundaryField()[patchi], patchi);
    }

//  Info << " finished " << endl;
    return tCv;
}


template<class MixtureType>
bool Foam::ePsiMixtureThermo<MixtureType>::read()
{
    if (eCombustionThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


template<class MixtureType>
void Foam::ePsiMixtureThermo<MixtureType>::eFromT(volScalarField& eToCompute, const volScalarField& Tgiven) const
{
      //notImplemented("ePsiMixtureThermo<MixtureType>::eFromT(volScalarField& Tgiven) const");
      //return volScalarField::null();
  if(Tgiven.size() == eToCompute.size())
  {
    forAll(Tgiven, celli)
    {
        eToCompute[celli] = this->cellMixture(celli).E(Tgiven[celli]);	// take the actual mixture composition, but at the temperature given by the user
    }
  }
  else
  {
            FatalErrorIn
            (
                "ePsiMixtureThermo<MixtureType>::eFromT "
                "(volScalarField& eToCompute, const volScalarField& Tgiven) const"
            )   << "size of energy field and temperature field provided are not the same"
                << abort(FatalError);
  }
  
}

template<class MixtureType>
void Foam::ePsiMixtureThermo<MixtureType>::hFromT(volScalarField& hToCompute, const volScalarField& Tgiven) const
{
  if(Tgiven.size() == hToCompute.size())
  {
    forAll(Tgiven, celli)
    {
        hToCompute[celli] = this->cellMixture(celli).H(Tgiven[celli]);	// take the actual mixture composition, but at the temperature given by the user
    }
  }
  else
  {
            FatalErrorIn
            (
                "ePsiMixtureThermo<MixtureType>::hFromT "
                "(volScalarField& hToCompute, const volScalarField& Tgiven) const"
            )   << "size of energy field and temperature field provided are not the same"
                << abort(FatalError);
  }
  
}


// ************************************************************************* //
