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
    namespace photoBio
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

Foam::photoBio::wideBandExtinction::wideBandExtinction
(
    const dictionary& dict
)
:
    extinctionModel(dict),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    absorption_(readBool(coeffsDict_.lookup("absorption"))),
    scattering_(readBool(coeffsDict_.lookup("scattering"))),
    nBands_(readLabel(coeffsDict_.lookup("nBands")))
{
    // Allocate absorption coefficient list (and initialize in case
    // absorption is turned off)
    aBand_.setSize(nBands_);
    forAll(aBand_, i)
    {
        aBand_[i] = 0.0;
    }
    
    // Allocate scattering coefficient list (and initialize in case
    // scattering is turned off)
    sBand_.setSize(nBands_);
    forAll(aBand_, i)
    {
        sBand_[i] = 0.0;
    }

    // Read the coefficients
    if (absorption_) coeffsDict_.lookup("absorptionCoeff") >> aBand_;
    if (scattering_) coeffsDict_.lookup("scatteringCoeff") >> sBand_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::photoBio::wideBandExtinction::~wideBandExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::photoBio::wideBandExtinction::a( const label bandI ) const
{        
	return  aBand_[bandI];
}


Foam::scalar Foam::photoBio::wideBandExtinction::s( const label bandI ) const
{        
	return  sBand_[bandI];
}


// ************************************************************************* //


