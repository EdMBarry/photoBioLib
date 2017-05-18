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

Class
    Foam::photoBio::extinctionModel

Description
    Model to supply absorption and scattering coefficients

\*---------------------------------------------------------------------------*/

#ifndef photoBioExtinctionModel_H
#define photoBioExtinctionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

/*---------------------------------------------------------------------------*\
                   Class absorptionEmissionModel Declaration
\*---------------------------------------------------------------------------*/

class extinctionModel
{

protected:

    // Protected data

        //- photoBio model dictionary
        const dictionary dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("extinctionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            extinctionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );

    // Constructors

        //- Construct from components
        extinctionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Selector
    static autoPtr<extinctionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~extinctionModel();

    // Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }

            //- Absorption coefficient
            virtual scalar a(const label bandI) const = 0;

            //- Scattering coefficient
            virtual scalar s(const label bandI) const = 0;
            
            //- Correct Extinction coefficients
            virtual void correct() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

