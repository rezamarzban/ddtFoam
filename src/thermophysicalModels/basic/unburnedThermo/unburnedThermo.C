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




#include "unburnedThermo.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedUnburnedEnthalpyFvPatchScalarField.H"
#include "gradientUnburnedEnthalpyFvPatchScalarField.H"
#include "mixedUnburnedEnthalpyFvPatchScalarField.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(unburnedThermo, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


Foam::wordList Foam::unburnedThermo::huBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf = Tu_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedUnburnedEnthalpyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientUnburnedEnthalpyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedUnburnedEnthalpyFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


void Foam::unburnedThermo::huBoundaryCorrection(volScalarField& hu)
{
    volScalarField::GeometricBoundaryField& hbf = hu.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientUnburnedEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientUnburnedEnthalpyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedUnburnedEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedUnburnedEnthalpyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Foam::unburnedThermo::unburnedThermo(const fvMesh& mesh, const volScalarField& fH)
Foam::unburnedThermo::unburnedThermo(const fvMesh& mesh, const volScalarField& p, const volScalarField& fH, const scalar& yO2inAir)
:
    IOdictionary
    (
        IOobject
        (
            "uthermoProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    yO2inAir_(yO2inAir),
    MH2_(2.01594),  // masses from atomicWeights.H
    MO2_(31.9988),
    MN2_(28.0134),
    Runiv_(8314.51),
    Tlow_(200.0),
    Thigh_(5000.0),
    Tcommon_(1000.0),
    maxIter_(50),
    Ttol_(0.05),
    Tmin_(270.0), // bounding
    Tmax_(2000.0), // bounding
    

    Tu_
    (
        IOobject
        (
            "Tu",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    hu_
    (
        IOobject
        (
            "hu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("hu_init",dimensionSet(0,2,-2,0,0),0.0),  // value is corrected in initialize/correct function
        this->huBoundaryTypes()
    ),
    p_(p),	// reference to the top-level pressure
    fH_(fH),	// reference to the top-level mixture fraction
    Ru_
    (
        IOobject
        (
            "Ru",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("Ru_init",dimensionSet(0,-2,2,-1,0),1e-20)  // value is corrected in initialize/correct function
    ),
    psiu_
    (
        IOobject
        (
            "psiu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("psiu_init",dimensionSet(0,-2,2,0,0),0.0)  // value is corrected in initialize/correct function
    ),
    rhou_
    (
        IOobject
        (
            "rhoUnburned",	// to distinguish it from rho*velocity
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("rhou_init",dimensionSet(1,-3,0,0,0),0.0)  // value is corrected in initialize/correct function
    )
    

{
    H2coeffs_[0]= 0.02991423E+02;
    H2coeffs_[1]= 0.07000644E-02;
    H2coeffs_[2]=-0.05633829E-06;
    H2coeffs_[3]=-0.09231578E-10;
    H2coeffs_[4]= 0.01582752E-13;
    H2coeffs_[5]=-0.08350340E+04;
    H2coeffs_[6]=-0.01355110E+02;

    H2coeffs_[7] = 0.03298124E+02;
    H2coeffs_[8] = 0.08249442E-02;
    H2coeffs_[9] =-0.08143015E-05;
    H2coeffs_[10]=-0.09475434E-09;
    H2coeffs_[11]= 0.04134872E-11;
    H2coeffs_[12]=-0.01012521E+05;
    H2coeffs_[13]=-0.03294094E+02;

    O2coeffs_[0]= 0.03697578E+02;
    O2coeffs_[1]= 0.06135197E-02;
    O2coeffs_[2]=-0.01258842E-05;
    O2coeffs_[3]= 0.01775281E-09;
    O2coeffs_[4]=-0.01136435E-13;
    O2coeffs_[5]=-0.01233930E+05;
    O2coeffs_[6]= 0.03189166E+02;

    O2coeffs_[7] = 0.03212936E+02;
    O2coeffs_[8] = 0.01127486E-01;
    O2coeffs_[9] =-0.05756150E-05;
    O2coeffs_[10]= 0.01313877E-07;
    O2coeffs_[11]=-0.08768554E-11;
    O2coeffs_[12]=-0.01005249E+05;
    O2coeffs_[13]= 0.06034738E+02;

    N2coeffs_[0]= 0.02926640E+02;
    N2coeffs_[1]= 0.01487977E-01;
    N2coeffs_[2]=-0.05684761E-05;
    N2coeffs_[3]= 0.01009704E-08;
    N2coeffs_[4]=-0.06753351E-13;
    N2coeffs_[5]=-0.09227977E+04;
    N2coeffs_[6]= 0.05980528E+02;

    N2coeffs_[7] = 0.03298677E+02;
    N2coeffs_[8] = 0.01408240E-01;
    N2coeffs_[9] =-0.03963222E-04;
    N2coeffs_[10]= 0.05641515E-07;
    N2coeffs_[11]=-0.02444855E-10;
    N2coeffs_[12]=-0.01020900E+05;
    N2coeffs_[13]= 0.03950372E+02;

    Info << "Creating unburnedThermo with Tu bounded to min. " << Tmin_ << " K" << endl;
    
    Tu_.max(Tmin_);	// bounding    
    Tu_.write();

    initialize();
    
    rhou_.write();
    //hu_.write();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::unburnedThermo::~unburnedThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




void Foam::unburnedThermo::initialize()
{
	bool udebug=true;


    forAll(hu_,i)
    {
	hu_[i] = hu(Tu_[i],fH_[i]);
    }
    hu_.correctBoundaryConditions();


    volScalarField Tu_before=Tu_;
   
    correct();
    
    volScalarField Tu_diff = mag(Tu_before - Tu_); // the difference before and after correction;
    
    if(max(Tu_diff).value() > 1.1*Ttol_)
	{
	  volScalarField("Tudiff",Tu_before-Tu_).write();
	  
            FatalErrorIn
            (
                "unburnedThermo::inititialize()"
            )   << "Calculating Tu -> hu -> Tu failed" << nl
		<< "maximum difference in Tu is " << max(Tu_diff).value() << " K, " << nl
                << abort(FatalError);
	}


    hu_.oldTime();
    psiu_.oldTime();


    if(udebug)
    {
	Info<< "unburnedThermo initialized:" << nl << "   hu  min...max   = " << min(hu_).value() << " ... " << max(hu_).value() << endl;
	Info<<                                        "   Tu  min...max   = " << min(Tu_).value() << " ... " << max(Tu_).value() << endl;
//	Info<<                                        "   psiu  min...max   = " << min(psiu_).value() << " ... " << max(psiu_).value() << endl;
	Info<<                                        "   rhou  min...max   = " << min(rhou_).value() << " ... " << max(rhou_).value() << endl;
//	Info<< 					      "   max. Tu difference after correction: " << max(Tu_orig).value() << endl;
	Info << endl;
    }
}





void Foam::unburnedThermo::correct()  // update Tu_ => Ru_, psiu_, rhou_
{
     const scalarField& fHCells = fH_.internalField();
     const scalarField& huCells = hu_.internalField();
     const scalarField& pCells = p_.internalField();     
     scalarField& TuCells = Tu_.internalField();
     scalarField& RuCells = Ru_.internalField();          
     scalarField& psiuCells = psiu_.internalField();
     scalarField& rhouCells = rhou_.internalField();
     
     //const scalar underrelax = 0.9; // underrelaxation <1; recommended due to nonlinear T-h curve of JANAF species
     const scalar underrelax = 0.99; // underrelaxation <1; necessary due to nonlinear T-h curve of JANAF species     

// iterative determination of temperature from polynom
  forAll(TuCells,celli)
  {
    scalar Test = TuCells[celli];
    scalar Tnew = TuCells[celli];
    int iter = 0;
    const scalar& fHi = fHCells[celli];

    do  
    {
        Test = Tnew;
        Tnew = Test - underrelax*(hu(Test,fHi) - huCells[celli])/Cpu(Test,fHi);
	
        if(Tnew<Tmin_)
	{
	    Tnew=Tmin_;
	    Test=Tmin_;
	    //Info << "Bounding Tu to " << Tmin << "K in cell " << celli << endl;
	}
        if(Tnew>Tmax_)
	{
	    Tnew=Tmax_;
	    Test=Tmax_;
	}
	
	

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "unburnedThermo::correct(const volScalarField& fH)"
            )   << "Maximum number of iterations exceeded" << nl
		<< "Test = " << Test << " K, " << nl
		<< "hu_ = " << huCells[celli] << " J/kg, h(Test,fH) = " << hu(Test,fHi) << " J/kg"
                << abort(FatalError);
        }


    } while (mag(Tnew - Test) > Ttol_);

    TuCells[celli]=Tnew;
    RuCells[celli] = Runiv_*(fHi/MH2_ + yO2inAir_*(1.0-fHi)/MO2_ + (1.0-yO2inAir_)*(1.0-fHi)/MN2_);
    psiuCells[celli] = 1.0/(RuCells[celli]*TuCells[celli]);
    rhouCells[celli] = pCells[celli]*psiuCells[celli];
  }
  
  forAll(Tu_.boundaryField(),patchi)
  {
         const fvPatchScalarField& pfH = fH_.boundaryField()[patchi];
         const fvPatchScalarField& phu = hu_.boundaryField()[patchi];	 
         const fvPatchScalarField& pp = p_.boundaryField()[patchi];	 	 
	 fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];
	 fvPatchScalarField& pRu = Ru_.boundaryField()[patchi];
	 fvPatchScalarField& ppsiu = psiu_.boundaryField()[patchi]; 
	 fvPatchScalarField& prhou = rhou_.boundaryField()[patchi]; 	 

         forAll(pTu, facei)
         {
	    scalar Test = pTu[facei];
	    scalar Tnew = pTu[facei];
	    int iter = 0;
	    const scalar& fHi = pfH[facei];

	    do  
	    {
		Test = Tnew;
		Tnew = Test - underrelax*(hu(Test,fHi) - phu[facei])/Cpu(Test,fHi);
	
		if(Tnew<Tmin_)
		{
		    //Info << "Bounding Tu to " << Tmin << "K in face " << facei << " on patch " << patchi << endl;		  
		    Tnew=Tmin_;
		    Test=Tmin_;
		}

		if (iter++ > maxIter_)
		{
		  FatalErrorIn
		  (
		    "unburnedThermo::correct(const volScalarField& fH)"
		  )   << "Maximum number of iterations exceeded" << nl
		      << "Test = " << Test << " K, " << nl
		      << "hu_ = " << phu[facei] << " J/kg, h(Test,fH) = " << hu(Test,fHi) << " J/kg"
		      << abort(FatalError);
		}
	    } while (mag(Tnew - Test) > Ttol_);

	    pTu[facei] = Tnew;
	    pRu[facei] = Runiv_*(fHi/MH2_ + yO2inAir_*(1.0-fHi)/MO2_ + (1.0-yO2inAir_)*(1.0-fHi)/MN2_);
	    ppsiu[facei] = 1.0/(pRu[facei]*pTu[facei]);
	    prhou[facei] = pp[facei]*ppsiu[facei];
	 }
  }


}

void Foam::unburnedThermo::correct(const volScalarField& T)  // update Tu_ => Ru_, psiu_, rhou_; Tu bounded by T
{
     const scalarField& fHCells = fH_.internalField();
     const scalarField& huCells = hu_.internalField();
     const scalarField& pCells = p_.internalField();     
     scalarField& TuCells = Tu_.internalField();
     scalarField& RuCells = Ru_.internalField();          
     scalarField& psiuCells = psiu_.internalField();
     scalarField& rhouCells = rhou_.internalField();
     const scalarField& TCells = T.internalField();
     
     const scalar underrelax = 1.0; // underrelaxation <1 recommended due to nonlinear T-h curve of JANAF species     

// iterative determination of temperature from polynom
  forAll(TuCells,celli)
  {
    scalar Tmax= TCells[celli];
    scalar Test = TuCells[celli];
    scalar Tnew = TuCells[celli];
    int iter = 0;
    const scalar& fHi = fHCells[celli];

    do  
    {
        Test = Tnew;
        Tnew = Test - underrelax*(hu(Test,fHi) - huCells[celli])/Cpu(Test,fHi);
	
        if(Tnew<Tmin_)
	{
	    Tnew=Tmin_;
	    Test=Tmin_;
	    //Info << "Bounding Tu to " << Tmin << "K in cell " << celli << endl;
	}
	if(Tnew>Tmax)
	{
	    Tnew=Tmax;
	    Test=Tmax;
	    //Info << "Bounding Tu to a maximum of " << Tmax << "K in cell " << celli << endl;
	}
	  
	

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "unburnedThermo::correct(const volScalarField& fH)"
            )   << "Maximum number of iterations exceeded" << nl
		<< "Test = " << Test << " K, " << nl
		<< "hu_ = " << huCells[celli] << " J/kg, h(Test,fH) = " << hu(Test,fHi) << " J/kg"
                << abort(FatalError);
        }


    } while (mag(Tnew - Test) > Ttol_);

    TuCells[celli]=Tnew;
    RuCells[celli] = Runiv_*(fHi/MH2_ + yO2inAir_*(1.0-fHi)/MO2_ + (1.0-yO2inAir_)*(1.0-fHi)/MN2_);
    psiuCells[celli] = 1.0/(RuCells[celli]*TuCells[celli]);
    rhouCells[celli] = pCells[celli]*psiuCells[celli];
  }
  
  forAll(Tu_.boundaryField(),patchi)
  {
         const fvPatchScalarField& pfH = fH_.boundaryField()[patchi];
         const fvPatchScalarField& phu = hu_.boundaryField()[patchi];	 
         const fvPatchScalarField& pp = p_.boundaryField()[patchi];	 	 
	 fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];
	 fvPatchScalarField& pRu = Ru_.boundaryField()[patchi];
	 fvPatchScalarField& ppsiu = psiu_.boundaryField()[patchi]; 
	 fvPatchScalarField& prhou = rhou_.boundaryField()[patchi]; 	 
	 const fvPatchScalarField& pT = T.boundaryField()[patchi];
	 
         forAll(pTu, facei)
         {
	    scalar Tmax = pT[facei];
	    scalar Test = pTu[facei];
	    scalar Tnew = pTu[facei];
	    int iter = 0;
	    const scalar& fHi = pfH[facei];

	    do  
	    {
		Test = Tnew;
		Tnew = Test - underrelax*(hu(Test,fHi) - phu[facei])/Cpu(Test,fHi);
	
		if(Tnew<Tmin_)
		{
		    Tnew=Tmin_;
		    Test=Tmin_;
		}
		if(Tnew>Tmax)
		{
		    Tnew=Tmax;
		    Test=Tmax;
		}

		if (iter++ > maxIter_)
		{
		  FatalErrorIn
		  (
		    "unburnedThermo::correct(const volScalarField& fH)"
		  )   << "Maximum number of iterations exceeded" << nl
		      << "Test = " << Test << " K, " << nl
		      << "hu_ = " << phu[facei] << " J/kg, h(Test,fH) = " << hu(Test,fHi) << " J/kg"
		      << abort(FatalError);
		}
	    } while (mag(Tnew - Test) > Ttol_);

	    pTu[facei] = Tnew;
	    pRu[facei] = Runiv_*(fHi/MH2_ + yO2inAir_*(1.0-fHi)/MO2_ + (1.0-yO2inAir_)*(1.0-fHi)/MN2_);
	    ppsiu[facei] = 1.0/(pRu[facei]*pTu[facei]);
	    prhou[facei] = pp[facei]*ppsiu[facei];
	 }
  }


}


Foam::scalar Foam::unburnedThermo::hu(const scalar& Tui, const scalar& fHi) const
{
	label move(0);	
	if(Tui<Tcommon_) move=7; else move=0;

	scalar hH2 = H2coeffs_[move+5];
	for(label j=0; j<5; j++)
	{
	    hH2 += H2coeffs_[move+j]/(j+1)*std::pow(Tui,j+1);
	}
	hH2 *= Runiv_/MH2_;

	scalar hO2 = O2coeffs_[move+5];
	for(label j=0; j<5; j++)
	{
	    hO2 += O2coeffs_[move+j]/(j+1)*std::pow(Tui,j+1);
	}
	hO2 *= Runiv_/MO2_;

	scalar hN2 = N2coeffs_[move+5];
	for(label j=0; j<5; j++)
	{
	    hN2 += N2coeffs_[move+j]/(j+1)*std::pow(Tui,j+1);
	}
	hN2 *= Runiv_/MN2_;

	/*
	Info << "   T = " << Tui << endl;
	Info << "   h: " << hH2 << " , " << hO2 << " , " << hN2 << endl;
	Info << "   insgesamt: " << fHi* hH2 + yO2inAir_*(1.0-fHi)*hO2 + (1.0-yO2inAir_)*(1.0-fHi)*hN2 << endl;
	*/

	return fHi* hH2 + yO2inAir_*(1.0-fHi)*hO2 + (1.0-yO2inAir_)*(1.0-fHi)*hN2;


}

Foam::tmp<Foam::scalarField> Foam::unburnedThermo::hu
(
    const scalarField& Tu,
    const scalarField& fH,
    const label patchi
) const
{
    tmp<scalarField> thu(new scalarField(Tu.size()));
    scalarField& hunburned = thu();

    forAll(Tu, facei)
    {
        //hu[facei] = this->patchFaceMixture(patchi, facei).H(Tu[facei]);
	hunburned[facei] = hu(Tu[facei],fH[facei]);
    }

    return thu;
}

Foam::tmp<Foam::scalarField> Foam::unburnedThermo::hu
(
    const scalarField& Tu,
    const scalarField& fH,
    const labelList& cells
) const
{
    tmp<scalarField> thu(new scalarField(Tu.size()));
    scalarField& hunburned = thu();

    forAll(Tu, celli)
    {
        hunburned[celli] = hu(Tu[celli],fH[celli]);
    }

    return thu;

}


Foam::scalar Foam::unburnedThermo::Cpu(const scalar& Tui, const scalar& fHi) const
{
//	return scalar(1000.0);

	label move(0);	
	if(Tui<Tcommon_) move=7; else move=0;

	scalar cpH2 = H2coeffs_[move];
	for(label j=1; j<5; j++)
	{
	    cpH2 += H2coeffs_[move+j]*std::pow(Tui,j);
	}
	cpH2 *= Runiv_/MH2_;

	scalar cpO2 = O2coeffs_[move];
	for(label j=1; j<5; j++)
	{
	    cpO2 += O2coeffs_[move+j]*std::pow(Tui,j);
	}
	cpO2 *= Runiv_/MO2_;

	scalar cpN2 = N2coeffs_[move];
	for(label j=1; j<5; j++)
	{
	    cpN2 += N2coeffs_[move+j]*std::pow(Tui,j);
	}
	cpN2 *= Runiv_/MN2_;

	//Info << " cp = " << fHi* cpH2 + yO2inAir_*(1.0-fHi)*cpO2 + (1.0-yO2inAir_)*(1.0-fHi)*cpN2 << endl;

	return fHi* cpH2 + yO2inAir_*(1.0-fHi)*cpO2 + (1.0-yO2inAir_)*(1.0-fHi)*cpN2;

}

Foam::tmp<Foam::scalarField> Foam::unburnedThermo::Cpu
(
    const scalarField& Tu,
    const scalarField& fH,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(Tu.size()));

    scalarField& cp = tCp();

    forAll(Tu, facei)
    {
        cp[facei] = Cpu(Tu[facei],fH[facei]);
    }

    return tCp;

}


Foam::tmp<Foam::volScalarField> Foam::unburnedThermo::Cpu() const // specific isobaric heat capacity of unburned mixture
{

    const fvMesh& mesh = Tu_.mesh();

    tmp<volScalarField> tCpu
    (
        new volScalarField
        (
            IOobject
            (
                "Cpu",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cpu = tCpu();

    scalarField& cpuCells = cpu.internalField();
    const scalarField& TuCells = Tu_.internalField();
    const scalarField& fHCells = fH_.internalField();    

    forAll(TuCells, celli)
    {
        //cpuCells[celli] = this->cellMixture(celli).Cpu(TuCells[celli],fHCells[celli]);
	cpuCells[celli] = Cpu(TuCells[celli],fHCells[celli]);
    }

    forAll(Tu_.boundaryField(), patchi)
    {
        cpu.boundaryField()[patchi] = Cpu(Tu_.boundaryField()[patchi], fH_.boundaryField()[patchi], patchi);
    }

    return tCpu;
    
}

void Foam::unburnedThermo::huFromT(volScalarField& huToCompute, const volScalarField& Tgiven) const
{
  if(Tgiven.size() == huToCompute.size())
  {
    forAll(Tgiven, celli)
    {
        huToCompute[celli] = hu(Tgiven[celli],fH_[celli]);	// take the actual mixture composition, but at the temperature given by the user
    }
  }
  else
  {
            FatalErrorIn
            (
                "unburnedThermo::huFromT "
                "(volScalarField& huToCompute, const volScalarField& Tgiven) const"
            )   << "size of energy field and temperature field provided are not the same"
                << abort(FatalError);
  }
  
}



/*
const Foam::volScalarField& Foam::unburnedThermo::p() const
{
    return p_;
}

const Foam::volScalarField& Foam::unburnedThermo::psiu() const
{
    return psiu_;
}


Foam::volScalarField& Foam::unburnedThermo::h()
{
    notImplemented("unburnedThermo::h()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::unburnedThermo::h() const
{
    notImplemented("unburnedThermo::h() const");
    return volScalarField::null();
}







Foam::volScalarField& Foam::unburnedThermo::hs()
{
    notImplemented("unburnedThermo::hs()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::unburnedThermo::hs() const
{
    notImplemented("unburnedThermo::hs() const");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::unburnedThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "unburnedThermo::hs"
        "(const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(NULL);
}


Foam::tmp<Foam::scalarField> Foam::unburnedThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "unburnedThermo::hs"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(NULL);
}


Foam::tmp<Foam::volScalarField> Foam::unburnedThermo::hc() const
{
    notImplemented("unburnedThermo::hc()");
    return volScalarField::null();
}


Foam::volScalarField& Foam::unburnedThermo::e()
{
    notImplemented("unburnedThermo::e()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::unburnedThermo::e() const
{
    notImplemented("unburnedThermo::e()");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::unburnedThermo::e
(
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "unburnedThermo::e"
        "(const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(NULL);
}


Foam::tmp<Foam::scalarField> Foam::unburnedThermo::e
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "unburnedThermo::e"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(NULL);
}


const Foam::volScalarField& Foam::unburnedThermo::T() const
{
    return T_;
}






Foam::tmp<Foam::scalarField> Foam::unburnedThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "unburnedThermo::Cv"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(NULL);
}


Foam::tmp<Foam::volScalarField> Foam::unburnedThermo::Cv() const
{
    notImplemented("unburnedThermo::Cv() const");
    return volScalarField::null();
}


const Foam::volScalarField& Foam::unburnedThermo::mu() const
{
    return mu_;
}


const Foam::volScalarField& Foam::unburnedThermo::alpha() const
{
    return alpha_;
}


bool Foam::unburnedThermo::read()
{
    return regIOobject::read();
}
*/


// ************************************************************************* //
