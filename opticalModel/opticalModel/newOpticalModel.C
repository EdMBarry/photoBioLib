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
#include "opticalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<opticalModel> opticalModel::New
(
     const volScalarField& T
)
{
    word opticalModelTypeName;

    // Note: no need to register/keep opticalProperties since models read
    // it themselves.
    {
        IOdictionary opticalPropertiesDict
        (
            IOobject
            (
                "opticalProperties",
                T.time().constant(),
                T.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        opticalPropertiesDict.lookup("opticalModel")
            >> opticalModelTypeName;
    }

    Info<< "Selecting opticalModel " << opticalModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(opticalModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "opticalModel::New(const volScalarField&)"
        )   << "Unknown opticalModel type " << opticalModelTypeName
            << nl << nl
            << "Valid opticalModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<opticalModel>(cstrIter()(T));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End optical
} // End namespace Foam

// ************************************************************************* //
