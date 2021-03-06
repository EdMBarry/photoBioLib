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
    Foam::photoBio::photoBioIntensityRay

Description
    Radiation intensity for a ray in a given direction

SourceFiles
    photoBioIntensityRay.C

\*---------------------------------------------------------------------------*/

#ifndef photoBioIntensityRay_H
#define photoBioIntensityRay_H

#include "extinctionModel.H"
#include "phaseFunctionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

// Forward declaration of classes
class photoBioDOM;

/*---------------------------------------------------------------------------*\
                    Class photoBioIntensityRay Declaration
\*---------------------------------------------------------------------------*/

class photoBioIntensityRay
{
public:

    static const word intensityPrefix;

private:

    // Private data

        //- Refence to the owner photoBioDOM object
        const photoBioDOM& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Ray intensity [W/m^2]
        autoPtr<volScalarField>  I_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Wavelength band index
        label iBand_;

        //- Angle index
        label iAngle_;

        //- Global ray id - incremented in constructor
        static label rayId;

    // Private member functions

        //- Disallow default bitwise copy construct
        photoBioIntensityRay(const photoBioIntensityRay&);

        //- Disallow default bitwise assignment
        void operator=(const photoBioIntensityRay&);

public:

    // Constructors

        //- Construct from components
        photoBioIntensityRay
        (
            const photoBioDOM& dom,
            const fvMesh& mesh,
            const label  iBand,
            const label  iAngle,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta

        );

    // Destructor
    ~photoBioIntensityRay();

    // Member functions

        // Edit

            //- Update photoBio intensity on i direction
            scalar correct();

            //- Update boundaries
            void updateBoundary();

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

        // Access

            //- Return intensity
             inline const volScalarField& I() const;

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

             //- Return the index of band
            inline  label iBand() const;

            //- Return the index of the angle
            inline label iAngle() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "photoBioIntensityRayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
