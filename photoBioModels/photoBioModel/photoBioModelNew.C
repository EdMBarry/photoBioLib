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
#include "photoBioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<photoBioModel> photoBioModel::New
(
     const volScalarField& Enu
)
{
    word photoBioModelTypeName;

    // Note: no need to register/keep photoBioProperties since models read
    // it themselves.
    {
        IOdictionary photoBioPropertiesDict
        (
            IOobject
            (
                "photoBioProperties",
                Enu.time().constant(),
                Enu.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        photoBioPropertiesDict.lookup("photoBioModel")
            >> photoBioModelTypeName;
    }

    Info<< "Selecting photoBioModel " << photoBioModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(photoBioModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "photoBioModel::New(const volScalarField&)"
        )   << "Unknown photoBioModel type " << photoBioModelTypeName
            << nl << nl
            << "Valid photoBioModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<photoBioModel>(cstrIter()(Enu));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End photoBio
} // End namespace Foam

// ************************************************************************* //
