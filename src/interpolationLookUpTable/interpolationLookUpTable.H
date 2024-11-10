/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::interpolationLookUpTable

Description
    A list of lists. Interpolates in all dimensions.
    The values must be monotonically increasing or decreasing in each dimension

Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element always returns the same value.
    - interpOutput_ contains now only the return values, not the values asked for
    - lookUp only for a scalar list (any dimension)

SourceFiles
    interpolationLookUpTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationLookUpTable_H
#define interpolationLookUpTable_H

#include "List.H"
#include "ListOps.H"
#include "scalarField.H"
#include "HashTable.H"
#include "IOdictionary.H"
#include "fvCFD.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class interpolationLookUpTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationLookUpTable
:
    public List<scalarField>
{
private:

    // Privsate data

        //- File name
        fileName fileName_;

        //- Table dimensions
        List<label> dim_;

        //- Min in each dimension
        List<scalar> min_;

        //- Maximum in each dimension
        List<scalar> max_;

	// coordinate is increasing (1) or decreasing (-1)
	List<label> increasing_;

        //- Dictionary entries
        List<dictionary> entries_;

        //- Output dictionaries
        List<dictionary> output_;

        //- Input indices from the look up table
        List<label> entryIndices_;

        //- Output Indices from the Look Up Table
        List<label> outputIndices_;

        //- Field names and indices
        HashTable<label> fieldIndices_;

        //- Output list containing input and interpolation values of outputs
        List<scalar> interpOutput_;

	//- The lower and upper values between the which the return values should remain
	// the bounded version "lookUpbounded" ensures that they actually remain in between
	List<scalar> lowerValues_;
	List<scalar> upperValues_;
	
	bool logarithmic_;



    // Private Member Functions

        //- Read the table of data from file
        void readTable(const word& instance, const fvMesh& mesh);
	//void readTable(const fileName& path, const fvMesh& mesh);
	void readTable(const IOdictionary& control);

        //- Dimension table from dictionaries input and output
        void dimensionTable();

        //- Find table index by scalarList
        label index(const List<scalar>&, bool) const;

        //- Check range of lookup value
        bool checkRange(const scalar, const label) const;

        //- Interpolate function return an scalar
        scalar interpolate
        (
            const label lo,
            const label hi,
            const scalar lookUpValue,
            const label ofield,
            const label interfield,
	    bool explain
        ) const;

        // Check list is monotonically increasing
        void check() const;

        // find hi index, interpolate and populate interpOutput_
        void findHi(const label lo, const List<scalar>& retvals, bool explain);



public:

    // Constructors

        //- Construct null
        interpolationLookUpTable();

        //- Construct given the name of the file containing the table of data
        interpolationLookUpTable
        (
            const fileName& fn,
            const word& instance,
            const fvMesh& mesh
        );

	interpolationLookUpTable
        (
            const IOdictionary& dict
        );
	

         //- Construct from dictionary
        interpolationLookUpTable(const dictionary& dict);

        //- Construct copy
        interpolationLookUpTable(const interpolationLookUpTable& interpTable);


    // Member Functions
	void info() const;

        //- Return true if the filed exists in the table
        bool found(const word& fieldName) const;

        //- Return the output list given an input scalar list
        const List<scalar>& lookUp(const List<scalar>&, bool);

        //- Return the output list given an input scalar list
        const List<scalar>& lookUpBounded(const List<scalar>&, bool);


        //- Write Look Up Table to filename.
        void write
        (
            Ostream& os,
            const fileName& fn,
            const word& instance,
            const fvMesh& mesh
        ) const;
	
	bool islogarithmic() const;


    // Access

        //- Return the index of a field by name
        inline label findFieldIndex(const word& fieldName) const;

        //- Return const access to the output dictionaries
        inline const List<dictionary>& output() const;

        //- Return const access tp the dictionary entries
        inline const List<dictionary>& entries() const;

        //- Return const access to the list of min dimensions
        inline const List<scalar>& min() const;

        //- Return const access to the list of dimensions
        inline const List<label>& dim() const;

        //- Return const access to the deltas in each dimension
        //inline const List<scalar>& delta() const;

        //- Return const access to the list of max dimensions
        inline const List<scalar>& max() const;

        //- Return const access to the information if the entries is increasing or decreasing (1/-1)
        inline const List<label>& increasing() const;


        //- Return const access to the table name
        inline word tableName() const;


     // Member Operators

        //- Return an element of constant List<scalar, Type>
        const scalarField& operator[](const label) const;

        //- Return an element of List<scalar, Type>
        scalarField& operator[](const label);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationLookUpTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template <class Type>
Foam::label Foam::interpolationLookUpTable<Type>::index
(
    const List<scalar>& values,
    bool explain
) const
{
    label totalIndex=0;

    forAll(dim_, i)
    {

        label dim = 1;
        for (int j = i + 1; j < dim_.size(); j++)
        {
            dim *= dim_[j] + 1;		// number of lines (entries) in the 1D table
        }
	
//	Info << "index(): i = " << i << ", dim = " << dim << endl;

        label tempIndex=dim_[i];

	for (label k=0; k<=dim_[i]; k++)
	{	// find the index of the highest entry that is still less or equal the requested value
	    if(values[i]*increasing_[i] <= List<scalarField>::operator[](i).operator[](k*dim)*increasing_[i])
	    {
		tempIndex = Foam::max(k-1,0);	// index in the current dimension i
		break;
	    }
	}


        totalIndex +=  dim  *Foam::min		// assemble all dimensions to find the index in the 1D table
            (
 //               Foam::max(label((values[i] - min_[i])/delta_[i]), 0),
 //               dim_[i]
               Foam::max(tempIndex, 0),
		dim_[i]  // 1 weniger als die tatsächliche Größe, damit noch Platz für hi bleibt
            );



	//Info << "   index(): totalIndex = " << totalIndex << endl;
    }


	if(explain)
	{
	  List<scalar> listfound(values.size(),1e10);
	  for(label j=0; j<entries_.size(); j++)
	  {
		listfound[j] = List<scalarField>::operator[](j).operator[](totalIndex);
	  }
	  Info <<  "index():   using " << listfound << " as basis" << " for request " << values <<  endl;
	}

    return totalIndex;

}

/*
template <class Type>
Foam::label Foam::interpolationLookUpTable<Type>::index
(
    const scalar value
) const
{
    label totalIndex = 0;
    //label i = 0;
    label dim = 1;

    label tempIndex=dim_[0]-1;

	for (label k=0; k<=dim_[0]; k++)
	{
	    Info << " k = " << k << ", entryValue = " << List<scalarField>::operator[](0).operator[](k*dim)*increasing_[0] << endl;
	    if(value*increasing_[0] <= List<scalarField>::operator[](0).operator[](k*dim)*increasing_[0])
	    {
		//Info << "   " << values[i]*increasing_[i] << " <= " << List<scalarField>::operator[](i).operator[](k*dim)*increasing_[i] << endl;
		tempIndex = Foam::max(k-1,0);
		break;
	    }
	}
	Info << "   index():  tempIndex = " << tempIndex << endl;


        totalIndex +=  dim  *Foam::min
        (
 //               Foam::max(label((values[i] - min_[i])/delta_[i]), 0),
 //               dim_[i]
               Foam::max(tempIndex, 0),
		dim_[0]  // immer 1 kleiner, damit noch Platz für hi bleibt
        );

    Info << "   index(): totalIndex = " << totalIndex << endl;

    return totalIndex;
}
*/

template<class Type>
bool Foam::interpolationLookUpTable<Type>::checkRange
(
    const scalar Value,
    const label interfield
) const
{
    if (Value*increasing_[interfield] >= min_[interfield]*increasing_[interfield] 
		&&  Value*increasing_[interfield] <= max_[interfield]*increasing_[interfield])
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<class Type>
void Foam::interpolationLookUpTable<Type>::info() const
{
        Info << "-----------" << endl;
	Info << " interpolationLookupTable: " << fileName_ << endl;
	Info << " dimensions: " << dim_ << endl;
	Info << " entries: " << entries_ << endl;
	Info << " output: " << output_ << endl;
	
	
	Info << " fieldIndices: " << fieldIndices_ << endl;
	Info << " min: " << min_ << endl;
	Info << " max: " << max_ << endl;
        Info << "-----------" << endl;
	
}

template<class Type>
Foam::scalar Foam::interpolationLookUpTable<Type>::interpolate
(
    const label lo,
    const label hi,
    const scalar Value,
    const label ofield,
    const label interfield,
    bool explain
) const
{
        if
        (
            List<scalarField>::operator[](interfield).operator[](hi)
         != List<scalarField>::operator[](interfield).operator[](lo)
        )
        {
            scalar output
            (
                List<scalarField>::operator[](ofield).operator[](lo)
              + (
                    List<scalarField>::operator[](ofield).operator[](hi)
                  - List<scalarField>::operator[](ofield).operator[](lo)
                )
               *(
                    Value
                  - List<scalarField>::operator[](interfield).operator[](lo)
                )
               /(
                    List<scalarField>::operator[](interfield).operator[](hi)
                  - List<scalarField>::operator[](interfield).operator[](lo)
                )
            );

	
	if(explain)
	{
	  scalar y0= 	List<scalarField>::operator[](ofield).operator[](lo);
	  scalar y1=    List<scalarField>::operator[](ofield).operator[](hi);
	  scalar x =	Value;
	  scalar x0=    List<scalarField>::operator[](interfield).operator[](lo);
	  scalar x1=    List<scalarField>::operator[](interfield).operator[](hi);
	  
	  Info << "    interpolating: " << "(" << y1 << " - " << y0 << ") * (" 
		<< x << " - " << x0 << ")/(" << x1 << " - " << x0 << ") = " << output-y0 << endl;
	}


            return output;
        }
        else
        {
	    Info << "interpolate: value[lo]=value[hi]" << endl;
            return List<scalarField>::operator[](ofield).operator[](lo);
        }
}


template<class Type>
void Foam::interpolationLookUpTable<Type>::dimensionTable()
{
    min_.setSize(entries_.size());
    dim_.setSize(entries_.size());
    //delta_.setSize(entries_.size());
    max_.setSize(entries_.size());
    increasing_.setSize(entries_.size());
    entryIndices_.setSize(entries_.size());
    outputIndices_.setSize(output_.size());
    label index = 0;
    label tableDim = 1;

    forAll(entries_,i)
    {
        dim_[i] = readLabel(entries_[i].lookup("N"));
        max_[i] = readScalar(entries_[i].lookup("max"));
        min_[i] = readScalar(entries_[i].lookup("min"));
	increasing_[i]= max_[i] > min_[i] ? 1 : (max_[i] < min_[i] ? -1 : 0);
	if(increasing_[i]==0)
        {	
		FatalErrorIn
        	(
            	"Foam::interpolationLookUpTable<Type>::dimensionTable()"
        	)   << "min and max values of coordinate " << i << " in " << tableName() << " are the same: "
		    << " no interpolation possible!" << nl << exit(FatalError);
	}
        //delta_[i] = (max_[i] - min_[i])/dim_[i];
        tableDim *= dim_[i] + 1;
        fieldIndices_.insert(entries_[i].lookup("name"), index);
        entryIndices_[i] = index;
        index++;

    }
    //Info << "dimensionTable(): increasing_ : " << increasing_ << endl;
    //Info << "dimensionTable(): entryIndices_ : " << entryIndices_ << endl;


    forAll(output_,i)
    {
        fieldIndices_.insert(output_[i].lookup("name"), index);
        outputIndices_[i] = index;
        index++;
    }

    List<scalarField>& internal = *this;

    internal.setSize(entries_.size() + output_.size());

    //interpOutput_.setSize(entries_.size() + output_.size());
    interpOutput_.setSize(output_.size());

    lowerValues_.setSize(output_.size());
    upperValues_.setSize(output_.size());

    forAll(internal, i)
    {
        internal[i].setSize(tableDim);
    }
}


template<class Type>
void Foam::interpolationLookUpTable<Type>::readTable
(
    const word& instance,
    const fvMesh& mesh
)
{
    IOdictionary control
    (
        IOobject
        (
            fileName_,
            instance,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    control.lookup("fields") >> entries_;
    control.lookup("output") >> output_;
    control.lookup("values") >> *this;

    dimensionTable();

    check();

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::readTable()"
        )   << "table is empty" << nl << exit(FatalError);
    }
}
/*
template<class Type>
void Foam::interpolationLookUpTable<Type>::readTable
(
    const fileName& path,
    const fvMesh& mesh
)
{
  Info << "reading table from complete path " << path << endl;
    IOdictionary control
    (
        IOobject
        (        
            path,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE

        )
    );

    control.lookup("fields") >> entries_;
    control.lookup("output") >> output_;
    control.lookup("values") >> *this;

    dimensionTable();

    check();

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::readTable()"
        )   << "table is empty" << nl << exit(FatalError);
    }
}
*/

template<class Type>
void Foam::interpolationLookUpTable<Type>::readTable
(
    const IOdictionary& dict
)
{
    dict.lookup("fields") >> entries_;
    dict.lookup("output") >> output_;
    dict.lookup("values") >> *this;
    
    //dict.lookup("logarithmic") >> logarithmic_;
    dictionary optionsDict = dict.subDict("options");
    logarithmic_ = optionsDict.lookupOrAddDefault<bool>("logarithmic",false);
    
    
    dimensionTable();

    check();

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::readTable()"
        )   << "table is empty" << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable()
:
    List<scalarField>(),
    fileName_("fileNameIsUndefined")
{}



template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
    const fileName& fn, const word& instance, const fvMesh& mesh
)
:
    List<scalarField>(),
    fileName_(fn),
    dim_(0),
    min_(0),
    //delta_(0.0),
    max_(0.0),
    increasing_(0),
    entries_(0),
    output_(0),
    entryIndices_(0),
    outputIndices_(0),
    interpOutput_(0),
    lowerValues_(0),
    upperValues_(0)
{
    readTable(instance, mesh);
}
/*
template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
    const fileName& path, const fvMesh& mesh
)
:
    List<scalarField>(),
    fileName_(path),
    dim_(0),
    min_(0),
    //delta_(0.0),
    max_(0.0),
    increasing_(0),
    entries_(0),
    output_(0),
    entryIndices_(0),
    outputIndices_(0),
    interpOutput_(0),
    lowerValues_(0),
    upperValues_(0)
{
    readTable(path, mesh);
}
*/

template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
    const IOdictionary& dict
)
:
    List<scalarField>(),
    fileName_(dict.name()),
    dim_(0),
    min_(0),
    //delta_(0.0),
    max_(0.0),
    increasing_(0),
    entries_(0),
    output_(0),
    entryIndices_(0),
    outputIndices_(0),
    interpOutput_(0),
    lowerValues_(0),
    upperValues_(0)
{
    Info << " Creating interpolationLookUpTable from dictionary" << endl;
    readTable(dict);
}



template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
     const interpolationLookUpTable& interpTable
)
:
    List<scalarField>(interpTable),
    fileName_(interpTable.fileName_),
    entryIndices_(interpTable.entryIndices_),
    outputIndices_(interpTable.outputIndices_),
    dim_(interpTable.dim_),
    min_(interpTable.min_),
    //delta_(interpTable.delta_),
    max_(interpTable.max_),
    increasing_(interpTable.increasing_),
    entries_(interpTable.entries_),
    output_(interpTable.entries_),
    interpOutput_(interpTable.interpOutput_),
    lowerValues_(interpTable.lowerValues_),
    upperValues_(interpTable.upperValues_)
{}


template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
    const dictionary& dict
)
:
    List<scalarField>(),
    fileName_(fileName(dict.lookup("fileName")).expand()),
    dim_(0),
    min_(0.0),
    //delta_(0.0),
    max_(0.0),
    increasing_(0),
    entries_(dict.lookup("fields")),
    output_(dict.lookup("output")),
    entryIndices_(0),
    outputIndices_(0),
    fieldIndices_(0),
    interpOutput_(0),
    lowerValues_(0),
    upperValues_(0)
{
    dimensionTable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationLookUpTable<Type>::check() const
{

//Info << "entries_ : " << entries_ << endl;

forAll(entries_,k)
{
    // check order in all coordinate dimensions
    scalar prevValue = List<scalarField>::operator[](k).operator[](0);
    label dim = 1;
    for (int j = k+1; j < dim_.size(); j++)
    {
        dim *= dim_[j] + 1;
    }

	//Info << "entry " << k << " , dim = " << dim << endl;

    for (label i = 1; i <= dim_[k]; i++)
    {
        label index = i*dim;
        const scalar currValue =
            List<scalarField>::operator[](k).operator[](index);

	//Info << "   index = " << index << ", currValue = " << currValue << endl;

	//Info << index << ":checkOrder() : " << prevValue << " -> " << currValue << endl;

        // avoid duplicate values (divide-by-zero error)
        if (increasing_[k]*currValue <= increasing_[k]*prevValue)
        {
            FatalErrorIn
            (
                "Foam::interpolationLookUpTable<Type>::checkOrder() const"
            )   << "out-of-order value: " << currValue 		
                << " at index " << index << nl 
		<< "(previous value:     " << prevValue 
                << " at index " << index-dim << ")"
		<< nl << exit(FatalError);
        }
        prevValue = currValue;
    }
}
}


template<class Type>
void Foam::interpolationLookUpTable<Type>::write
(
    Ostream& os,
    const fileName& fn,
    const word& instance,
    const fvMesh& mesh
) const
{
    IOdictionary control
    (
        IOobject
        (
            fn,
            instance,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    control.writeHeader(os);

    os.writeKeyword("fields");
    os << entries_ << token::END_STATEMENT << nl;

    os.writeKeyword("output");
    os << output_ << token::END_STATEMENT << nl;

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationTable<Type>::write()"
        )   << "table is empty" << nl << exit(FatalError);
    }
    os.writeKeyword("values");
    os << *this << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::scalarField&
Foam::interpolationLookUpTable<Type>::operator[](const label i)
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::operator[](const label)"
        )   << "table has (" << n << ") columns" << nl << exit(FatalError);
    }
    else if (ii < 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::operator[](const label)"
        )   << "index (" << ii << ") underflow" << nl << exit(FatalError);
    }
    else if (ii > n)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::operator[](const label)"
        )   << "index (" << ii << ") overflow" << nl << exit(FatalError);
    }

    return List<scalarField>::operator[](ii);
}


template<class Type>
const Foam::scalarField&
Foam::interpolationLookUpTable<Type>::operator[](const label i) const
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::operator[]"
            "(const label) const"
        )   << "table has (" << n << ") columns" << nl << exit(FatalError);
    }
    else if (ii < 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::operator[]"
            "(const label) const"
        )   << "index (" << ii << ") underflow" << nl << exit(FatalError);
    }

    else if (ii > n)
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::operator[]"
            "(const label) const"
        )   << "index (" << ii << ") overflow" << nl
            << exit(FatalError);
    }

    return List<scalarField>::operator[](ii);
}


template<class Type>
bool Foam::interpolationLookUpTable<Type>::found(const word& fieldName) const
{
    return fieldIndices_.found(fieldName);
}




template<class Type>
const Foam::scalarList&
Foam::interpolationLookUpTable<Type>::lookUpBounded(const List<scalar>& retvals, bool explain)
{
	List<scalar> unboundedOutput = lookUp(retvals,explain);
    if(explain)
	{	
		Info << " lookUpBounded: " << nl << lowerValues_ << nl << upperValues_ << nl << interpOutput_ << endl;
	}
    forAll(interpOutput_,j)
    {
	interpOutput_[j] = upperValues_[j] > lowerValues_[j] ? Foam::max(Foam::min(interpOutput_[j],upperValues_[j]),lowerValues_[j]) : Foam::min(Foam::max(interpOutput_[j],upperValues_[j]),lowerValues_[j]);
    }
    if(explain)
	{	
		Info << " lookUpBounded: " << interpOutput_ << endl;
	}

    return interpOutput_;
}

template<class Type>
const Foam::scalarList&
Foam::interpolationLookUpTable<Type>::lookUp(const List<scalar>& retvals, bool explain)
{
    if(retvals.size() != entries_.size())
    {
        FatalErrorIn
        (
            "Foam::interpolationLookUpTable<Type>::lookUp(const List<scalar>& retvals)" 
        )   << " the table has "<< entries_.size() << " dimension, but retvals has length " << retvals.size()
            << exit(FatalError);
    }

    label lo = index(retvals,explain);	// find the line in the table with the values "below"
    findHi(lo, retvals, explain);	// find the line in the table with the values "above" and interpolate
    return interpOutput_;
}



template<class Type>
void Foam::interpolationLookUpTable<Type>::findHi
(
    const label lo,
    const List<scalar>& retvals,
    bool explain
)
{
    List<label> hi(entryIndices_.size(),-1); //Josef
    //List<label> hi(output_.size(),-1);
    //List<scalar> baseValue(output_.size(),0.0);
    label maxhi = lo;

    forAll(outputIndices_,j)
    {
	//baseValue[j] = List<scalarField>::operator[](outputIndices_[j]).operator[](lo);
	lowerValues_[j] = List<scalarField>::operator[](outputIndices_[j]).operator[](lo);
    }
	if(explain)
	{    
		//Info << "findHi():  baseValue = " << baseValue << endl;
		Info << "findHi():  baseValue = " << lowerValues_ << endl;
	}

    forAll(outputIndices_,j)
    {
        scalar tmp = 0;
        label ofield = outputIndices_[j];

	maxhi=lo; //Josef
	//Info << "findHi(): outputIndices j = " << j << ", baseValue = " << baseValue << " (ofield = " << ofield << ", lo = " << lo << ")"<<  endl;

        forAll(entryIndices_,i)
        {
            if (checkRange(retvals[i], entryIndices_[i])) // if the variable requested is out of the table bounds, do nothing - tmp gets the lowerValues
            {

	        label dim = 1;
	        for (int k = i + 1; k < dim_.size(); k++)
	        {
        	    dim *= dim_[k] + 1;
        	}
		//Info << " entryIndices_[" << i << "]" << endl;
		//Info << "    dim = " << dim << endl;

		hi[i] = lo+dim;

		maxhi +=dim; //zur Info

		//Info << "    lo = " << lo << ", hi[" << i << "] = " << hi[i] << endl;

		//Info << "findHi(): entryIndex i = " << i << endl;
		scalar dtmp = interpolate(lo, hi[i], retvals[i], ofield, entryIndices_[i],explain)
                     - lowerValues_[j];
		tmp += dtmp;

		//Info << "  dtmp = " << dtmp<< endl;
            }
            //interpOutput_[entryIndices_[i]] = retvals[i];
        }

        tmp += lowerValues_[j];
	//Info << " tmp = " << tmp << endl;

        //interpOutput_[outputIndices_[j]] = tmp;
	interpOutput_[j] = tmp;
    }

    forAll(outputIndices_,j)
    {
	upperValues_[j] = List<scalarField>::operator[](outputIndices_[j]).operator[](maxhi);
    }


if(explain)
{
    List<scalar> upperCoordinates(entries_.size(),0.0);
    //List<scalar> upperValue(output_.size(),0.0);

    forAll(entries_,j)
    {
	upperCoordinates[j] = List<scalarField>::operator[](j).operator[](maxhi);
    }
    Info << "findHi(): using " << upperCoordinates << " as upper bound" << endl;

    Info << "findHi(): upperValue = " << upperValues_ << endl;
}


}

template<class Type>
bool Foam::interpolationLookUpTable<Type>::islogarithmic() const
{
    return logarithmic_;
}




// ************************************************************************* //


//#ifdef NoRepository
//#   include "interpolationLookUpTable.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
