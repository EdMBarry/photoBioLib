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

#include "wideBandExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(wideBandExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            wideBandExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::optical::wideBandExtinction::wideBandExtinction
(
    const dictionary& dict
 //   , const fvMesh& mesh
)
:
 //   extinctionModel(dict, mesh),
      extinctionModel(dict),
    coeffsDict_((dict.subDict(typeName + "Coeffs")))
{

    const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    
    functionDicts.lookup("absorption") >> absorption_;
    
    functionDicts.lookup("scattering") >>  scattering_;
    
    functionDicts.lookup("nBand") >> nBands_;
    
    functionDicts.lookup("bioDensity") >> bioDensity_;       
    
    aBand_.setSize(nBands_);  forAll(aBand_,i)  aBand_[i] = 0.0;
    sBand_.setSize(nBands_);  forAll(aBand_,i)  sBand_[i] = 0.0;
    
    if(absorption_)
    {
		functionDicts.lookup("absorptionCoeff") >> aBand_;
		forAll(aBand_,i)  aBand_[i] = aBand_[i]*bioDensity_;
	}

    if(scattering_)
    {
		functionDicts.lookup("scatteringCoeff") >> sBand_;
		forAll(sBand_,i)  sBand_[i] = sBand_[i]*bioDensity_;
	}
 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::wideBandExtinction::~wideBandExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
      
Foam::scalar Foam::optical::wideBandExtinction::a( const label bandI ) const
{        
	return  aBand_[bandI];
}
         
Foam::scalar Foam::optical::wideBandExtinction::s( const label bandI ) const
{        
	return  sBand_[bandI];
}



/*
 * 
 * Foam::scalar  Foam::optical::wideBandExtinction::k( const label bandI ) const
{
	return  kBand_[bandI];
}
           
           
Foam::tmp<Foam::volScalarField>
Foam::optical::wideBandExtinction::addIntensity
(
    const label i,
    const volScalarField& ILambda
) const
{
    return ILambda*(iBands_[i][1] - iBands_[i][0])/totalWaveLength_;
}



void Foam::optical::wideBandExtinction::correct
(
    volScalarField& k,
    PtrList<volScalarField>& kLambda
) const
{
    k = dimensionedScalar("zero", dimless/dimLength, 0.0);

    for (label j=0; j<nBands_; j++)
    {
        Info<< "Calculating extinction in band: " << j << endl;
        
        kLambda[j].internalField() = this->k(j);
        
        k.internalField() +=
            kLambda[j].internalField()
           *(iBands_[j][1] - iBands_[j][0])
           /totalWaveLength_;
    }

}

*/



