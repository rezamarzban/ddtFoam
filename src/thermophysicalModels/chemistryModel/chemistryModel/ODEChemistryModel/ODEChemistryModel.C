/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "ODEChemistryModel.H"
#include "chemistrySolver.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ODEChemistryModel<CompType, ThermoType>::ODEChemistryModel
(
    const fvMesh& mesh,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    CompType(mesh, thermoTypeName),

    ODE(),

    Y_(this->thermo().composition().Y()),

    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),

    solver_
    (
        chemistrySolver<CompType, ThermoType>::New
        (
            *this,
            compTypeName,
            thermoTypeName
        )
    ),

    RR_(nSpecie_)
{
    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );
    }

    Info<< "ODEChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ODEChemistryModel<CompType, ThermoType>::~ODEChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalarField Foam::ODEChemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalarField om(nEqns(), 0.0);

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return om;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c2(nSpecie_, 0.0);
    for (label i=0; i<nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(T, p, c2);
    scalar kr = R.kr(kf, T, p, c2);

    pf = 1.0;
    pr = 1.0;

    label Nl = R.lhs().size();
    label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        scalar exp = R.rhs()[srRef].exponent;
        if (exp<1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    scalar T = c[nSpecie_];
    scalar p = c[nSpecie_ + 1];

    dcdt = omega(c, T, p);

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        scalar W = specieThermo_[i].W();
        cSum += c[i];
        rho += W*c[i];
    }
    scalar mw = rho/cSum;
    scalar cp = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        scalar cpi = specieThermo_[i].cp(T);
        scalar Xi = c[i]/rho;
        cp += Xi*cpi;
    }
    cp /= mw;

    scalar dT = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        scalar hi = specieThermo_[i].h(T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    // limit the time-derivative, this is more stable for the ODE
    // solver when calculating the allowed time step
    scalar dtMag = min(500.0, mag(dT));
    dcdt[nSpecie_] = -dT*dtMag/(mag(dT) + 1.0e-10);

    // dp/dt = ...
    dcdt[nSpecie_+1] = 0.0;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    scalar T = c[nSpecie_];
    scalar p = c[nSpecie_ + 1];

    scalarField c2(nSpecie_, 0.0);
    for (label i=0; i<nSpecie_; i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    // length of the first argument must be nSpecie()
    dcdt = omega(c2, T, p);

    for (label ri=0; ri<reactions_.size(); ri++)
    {
        const Reaction<ThermoType>& R = reactions_[ri];

        scalar kf0 = R.kf(T, p, c2);
        scalar kr0 = R.kr(T, p, c2);

        forAll(R.lhs(), j)
        {
            label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c2[si]>SMALL)
                        {
                            kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c2[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c2[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar er = R.rhs()[i].exponent;
                if (i==j)
                {
                    if (er<1.0)
                    {
                        if (c2[si]>SMALL)
                        {
                            kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c2[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= pow(c2[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] -= sr*kr;
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;
    scalarField dcdT0 = omega(c2, T - delta, p);
    scalarField dcdT1 = omega(c2, T + delta, p);

    for (label i=0; i<nEqns(); i++)
    {
        dfdc[i][nSpecie()] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }

}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::tc() const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc();

    label nReaction = reactions_.size();

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];
            scalarField c(nSpecie_);
            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                cSum += c[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];

                omega
                (
                    R, c, Ti, pi, pf, cf, lRef, pr, cr, rRef
                );

                forAll(R.rhs(), s)
                {
                    scalar sr = R.rhs()[s].stoichCoeff;
                    tc[celli] += sr*pf*cf;
                }
            }
            tc[celli] = nReaction*cSum/tc[celli];
        }
    }


    ttc().correctBoundaryConditions();

    return ttc;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                scalar hi = specieThermo_[i].Hc();
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }

    return tSh;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class ThermoType>
Foam::label Foam::ODEChemistryModel<CompType, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::calculate()
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    for (label i=0; i<nSpecie_; i++)
    {
        RR_[i].setSize(rho.size());
    }

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0.0;
            }

            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];

            scalarField c(nSpecie_);
            scalarField dcdt(nEqns(), 0.0);

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
            }

            dcdt = omega(c, Ti, pi);

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
            }
        }
    }
}


template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::solve
(
    const scalar t0,
    const scalar deltaT
)
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    for (label i=0; i<nSpecie_; i++)
    {
        RR_[i].setSize(rho.size());
    }

    if (!this->chemistry_)
    {
        return GREAT;
    }
    
    if(deltaT<=SMALL)
    {
	  for (label i=0; i<nSpecie_; i++)
	  {
	    //Info << "max. RR[" << i << "]  vorher: " << max(RR_[i]) << endl;
	     RR_[i] *= 0.0;
	    //Info << "max. RR[" << i << "] nachher: " << max(RR_[i]) << endl;	     
	  }
        return GREAT;
    }


    scalar deltaTMin = GREAT;

    tmp<volScalarField> thc = this->thermo().hc();
    const scalarField& hc = thc();

    forAll(rho, celli)
    {
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = 0.0;
        }

        scalar rhoi = rho[celli];
        scalar Ti = this->thermo().T()[celli];
        scalar hi = this->thermo().hs()[celli] + hc[celli];
        scalar pi = this->thermo().p()[celli];

        scalarField c(nSpecie_);
        scalarField c0(nSpecie_);
        scalarField dc(nSpecie_, 0.0);

//scalar Yall(0.0);	
        for (label i=0; i<nSpecie_; i++)
        {
            c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
//Yall += Y_[i][celli];
        }
        c0 = c;
	

//Info << "Yall = " << Yall << endl;

        scalar t = t0;
        scalar tauC = this->deltaTChem_[celli];
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        // calculate the chemical source terms
        scalar cTot = sum(c);
	
        while (timeLeft > SMALL)
        {
            tauC = solver().solve(c, Ti, pi, t, dt);
            t += dt;
            // update the temperature
            cTot = sum(c);
            ThermoType mixture(0.0*specieThermo_[0]);
            for (label i=0; i<nSpecie_; i++)
            {
                mixture += (c[i]/cTot)*specieThermo_[i];
            }
            Ti = mixture.TH(hi, Ti);

            timeLeft -= dt;
            this->deltaTChem_[celli] = tauC;
            dt = min(timeLeft, tauC);
            dt = max(dt, SMALL);
        }
        deltaTMin = min(tauC, deltaTMin);

        dc = c - c0;

	/*
	// computation of WTot unnecessary?
        scalar WTot = 0.0;
        for (label i=0; i<nSpecie_; i++)
        {
            WTot += c[i]*specieThermo_[i].W();
        }
        WTot /= cTot;
	*/
	
	scalar RRges(0.0);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dc[i]*specieThermo_[i].W()/deltaT;
	    RRges += RR_[i][celli];
        }
if(abs(RRges)>100.0) // arbitrary limit
	{
	    Info << "ODEChemistryModel: The reaction rates don't cancel out in cell " << celli << ": " << RRges << " kg/m^3/s." << endl;
	    for (label i=0; i<nSpecie_; i++)
	    {
		Info << "  " << i << ": \t" << RR_[i][celli] << endl;
		//RR_[i][celli]=0.0;
	    }
	    Info << "  T:" << this->thermo().T()[celli] << " -> \t" << Ti << endl;
	}
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}

template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::solveIsobar
(
    const scalar t0,
    const scalar deltaT,
    UPtrList<volScalarField> RRates, 	// where to write the results
    const volScalarField& temperature,	// which temperature to use for evaluating the source terms
    const volScalarField& enthalpy	// the enthalpy which corresponds to the temperature -> does not change during reaction
)
{
  Info << "Using my own ODEChemistryModel::solveIsobar" << endl;
  
    if (!this->chemistry_)
    {
        return GREAT;
    }

    scalar deltaTMin = GREAT;

    forAll(temperature, celli)
    {
        for (label i=0; i<nSpecie_; i++)
        {
            RRates[i][celli] = 0.0;
        }
        
        scalar invWTot = 0.0;
        for (label i=0; i<nSpecie_; i++)
        {
            invWTot += Y_[i][celli]/specieThermo_[i].W();
        }
        scalar WTot = 1.0/invWTot;
	

        scalar Ti = temperature[celli]; 
	scalar hi = enthalpy[celli];		// remains constant (isobaric combustion)
        scalar pi = this->thermo().p()[celli]; 	// the same pressure in both parts of MinMax; remains constant (isobaric combustion)

        scalarField c(nSpecie_);
        scalarField c0(nSpecie_);
        scalarField dc(nSpecie_, 0.0);
	
	scalar rhoi = pi*WTot/(specieThermo_[0].RR*Ti);

        for (label i=0; i<nSpecie_; i++)
        {
	    c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
        }
        c0 = c;

        scalar t = t0;
        scalar tauC = this->deltaTChem_[celli];	// vorläufig // should be a different time scale on each side
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        scalar cTot = sum(c);
	
        while (timeLeft > SMALL)
        {
            tauC = solver().solve(c, Ti, pi, t, dt);
            t += dt;
            cTot = sum(c);

            ThermoType mixture(0.0*specieThermo_[0]);
            for (label i=0; i<nSpecie_; i++)
            {
                mixture += (c[i]/cTot)*specieThermo_[i];
            }
	    Ti = mixture.TH(hi, Ti);	// isobaric: evaluate new temperature for the enthalpy which remains constant

            timeLeft -= dt;
            this->deltaTChem_[celli] = tauC;
            dt = min(timeLeft, tauC);
            dt = max(dt, SMALL);
        }
        deltaTMin = min(tauC, deltaTMin);

        dc = c - c0;

	scalar RRges(0.0);

        for (label i=0; i<nSpecie_; i++)
        {
            RRates[i][celli] = dc[i]*specieThermo_[i].W()/deltaT;
	    RRges += RRates[i][celli];
        }
        
	if(abs(RRges)>100.0) // arbitrary limit
	{
	    Info << "ODEChemistryModel: The reaction rates don't cancel out in cell " << celli << ": " << RRges << " kg/m^3/s." << endl;
	    for (label i=0; i<nSpecie_; i++)
	    {
		Info << "  " << i << ": \t" << RR_[i][celli] << endl;
		//RR_[i][celli]=0.0;
	    }
	    Info << "  T:" << this->thermo().T()[celli] << " -> \t" << Ti << endl;
	}
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}



template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::solveIsochor // for MinMax model the respective temperature field has to be passed as well
(
    const scalar t0,
    const scalar deltaT,
    UPtrList<volScalarField> RRates, 	// where to write the results
    const volScalarField& temperature,	// which temperature to use for evaluating the source terms
    const volScalarField& energy	// the energy which corresponds to the temperature -> does not change during reaction
)
{
  Info << "Using my own ODEChemistryModel::solveIsochor" << endl;
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    if (!this->chemistry_)
    {
        return GREAT;
    }

    scalar deltaTMin = GREAT;

    //tmp<volScalarField> thc = this->thermo().hc(); 	
    //const scalarField& hc = thc();			

    forAll(rho, celli)
    {
        for (label i=0; i<nSpecie_; i++)
        {
            RRates[i][celli] = 0.0;
        }
        
        scalar invWTot = 0.0;
        for (label i=0; i<nSpecie_; i++)
        {
            invWTot += Y_[i][celli]/specieThermo_[i].W();
        }
        scalar WTot = 1.0/invWTot;
	

        //scalar rhoi = rho[celli];
        scalar Ti = temperature[celli]; 
	//Info << "solveIsochor: Ti = " << Ti ;
        //scalar hi = this->thermo().hs()[celli] + hc[celli];
	//scalar hi = this->thermo().h()[celli];
	//scalar hi = this->thermo().h(temperature[celli],celli);
	scalar ei = energy[celli];	// remains constant (isochoric combustion)
        scalar pi = this->thermo().p()[celli]; // the same pressure in both parts of MinMax
	//Info << "solveIsochor: pi = " << pi ;	

        scalarField c(nSpecie_);
        scalarField c0(nSpecie_);
        scalarField dc(nSpecie_, 0.0);
	
	scalar rhoi = pi*WTot/(specieThermo_[0].RR*Ti);

        for (label i=0; i<nSpecie_; i++)
        {
	    c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
	    //c[i] = pi*WTot/(specieThermo_[0].RR*Ti)*Y_[i][celli]/specieThermo_[i].W();
        }
        c0 = c;

        scalar t = t0;
        scalar tauC = this->deltaTChem_[celli];	// vorläufig // should be a different time scale on each side
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        // calculate the chemical source terms
        scalar cTot = sum(c);
	
/*	
	scalar WTot = 0.0; // specific molecular weight of the mixture
	    for (label i=0; i<nSpecie_; i++)	// check if pi fulfils ideal gas law
	    {
	      WTot += c[i]*specieThermo_[i].W();
	    }
	    WTot /= cTot;
	    //Info << "Rgas: " << specieThermo_[0].RR << endl;
	    //pi = rhoi*specieThermo_[0].RR/WTot*Ti;		// ideal gas law
*/
	    pi = cTot*specieThermo_[0].RR*Ti;			// ideal gas law
	 //   Info << " vs " << pi << endl;

	
	
        while (timeLeft > SMALL)
        {
            tauC = solver().solve(c, Ti, pi, t, dt);
            t += dt;
            cTot = sum(c);

            ThermoType mixture(0.0*specieThermo_[0]);
            for (label i=0; i<nSpecie_; i++)
            {
                mixture += (c[i]/cTot)*specieThermo_[i];
            }
            //Ti = mixture.TH(hi, Ti);
	    Ti = mixture.TE(ei, Ti);	// isochoric: evaluate new temperature for the internal energy which remains constant

	    // isochoric combustion requires update of pressure pi:
	    /*
	    WTot = 0.0;
	    for (label i=0; i<nSpecie_; i++)
	    {
	      WTot += c[i]*specieThermo_[i].W();
	    }
	    WTot /= cTot;
	    //Info << "Rgas: " << specieThermo_[0].RR << endl;
	    pi = rhoi*specieThermo_[0].RR/WTot*Ti;		// ideal gas law
	    */
	    pi = cTot*specieThermo_[0].RR*Ti;
	    //Info << "   " << pi << endl;		// ideal gas law

            timeLeft -= dt;
            this->deltaTChem_[celli] = tauC;
            dt = min(timeLeft, tauC);
            dt = max(dt, SMALL);
	    
        }
        deltaTMin = min(tauC, deltaTMin);

        dc = c - c0;
//Info << " -> " << Ti << endl;	

	scalar RRges(0.0);

        for (label i=0; i<nSpecie_; i++)
        {
            RRates[i][celli] = dc[i]*specieThermo_[i].W()/deltaT;
	    RRges += RRates[i][celli];
        }
if(abs(RRges)>100.0) // arbitrary limit
	{
	    Info << "ODEChemistryModel: The reaction rates don't cancel out in cell " << celli << ": " << RRges << " kg/m^3/s." << endl;
	    for (label i=0; i<nSpecie_; i++)
	    {
		Info << "  " << i << ": \t" << RR_[i][celli] << endl;
		//RR_[i][celli]=0.0;
	    }
	    Info << "  T:" << this->thermo().T()[celli] << " -> \t" << Ti << endl;
	}
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}

template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::updateMinMax // for MinMax model: add the reaction rates
(
    const volScalarField& alphaHI,
    UPtrList<volScalarField> RRatesHI,
    UPtrList<volScalarField> RRatesLO
)
{
  
	//alphaHI.min(1.0);
        for (label i=0; i<nSpecie_; i++)
        {
	  RR_[i] = alphaHI*RRatesHI[i] + (scalar(1.0)-alphaHI)*RRatesLO[i];
	}
	
}
// ************************************************************************* //
