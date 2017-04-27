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

#include "opticalModel.H"
#include "extinctionModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(opticalModel, 0);
        defineRunTimeSelectionTable(opticalModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::opticalModel::opticalModel(const volScalarField& intensity)
:
    IOdictionary
    (
        IOobject
        (
            "opticalProperties",
            intensity.time().constant(),
            intensity.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(intensity.mesh()),
    time_(intensity.time()),
    optical_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    extinction_(NULL)     
{}


Foam::optical::opticalModel::opticalModel
(
    const word& type,
    const volScalarField& intensity
)
:
    IOdictionary
    (
        IOobject
        (
            "opticalProperties",
            intensity.time().constant(),
            intensity.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(intensity.mesh()),
    time_(intensity.time()),
    optical_(lookup("optical")),
    coeffs_(subDict(type + "Coeffs")),
    solverFreq_(readLabel(lookup("solverFreq"))),
    extinction_(extinctionModel::New(*this))    
{
    solverFreq_ = max(1, solverFreq_);
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::optical::opticalModel::~opticalModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::optical::opticalModel::read()
{
    if (regIOobject::read())
    {
        lookup("optical") >> optical_;
        coeffs_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::optical::opticalModel::correct()
{
   if (!optical_)
    {
        return;
    }

    if (time_.timeIndex() % solverFreq_ == 0)
    {
        calculate();
    }
}




// ************************************************************************* //
