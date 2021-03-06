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

#include "error.H"
#include "phaseFunctionModel.H"
#include "photoBioDOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(phaseFunctionModel, 0);
        defineRunTimeSelectionTable(phaseFunctionModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::phaseFunctionModel::phaseFunctionModel
(
    const photoBioDOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    dom_(dom),
    dict_(dict),
    nDim_(nDim)
{}

// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::photoBio::phaseFunctionModel::~phaseFunctionModel()
{}


// ************************************************************************* //


bool Foam::photoBio::phaseFunctionModel::inScatter() const
{
    return false;
}


Foam::scalar Foam::photoBio::phaseFunctionModel::correct
(
    const label rayI,
    const label rayJ,
    const label iBand
)   const
{
    return 0.0;
}

