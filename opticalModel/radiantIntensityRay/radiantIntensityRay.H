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
    Foam::radiation::radiantIntensityRay

Description
    Radiation intensity for a ray in a given direction

SourceFiles
    radiantIntensityRay.C

\*---------------------------------------------------------------------------*/

#ifndef radiantIntensityRay_H
#define radiantIntensityRay_H

#include "absorptionEmissionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

// Forward declaration of classes
class opticalDOM;

/*---------------------------------------------------------------------------*\
                    Class radiantIntensityRay Declaration
\*---------------------------------------------------------------------------*/

class radiantIntensityRay
{
public:

    static const word intensityPrefix;


private:

    // Private data

        //- Refence to the owner opticalDOM object
        const opticalDOM& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Absorption/emission model
        const absorptionEmissionModel& absorptionEmission_;


        //- Total optical intensity / [W/m2]
        volScalarField I_;

        //- Total optical heat flux on boundary
   //     volScalarField Qr_;

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

        //- Number of wavelengths/bands
        label nLambda_;

        //- List of pointers to optical intensity fields for given wavelengths
        PtrList<volScalarField> ILambda_;

        //- Global ray id - incremented in constructor
        static label rayId;


    // Private member functions

        //- Disallow default bitwise copy construct
        radiantIntensityRay(const radiantIntensityRay&);

        //- Disallow default bitwise assignment
        void operator=(const radiantIntensityRay&);


public:

    // Constructors

        //- Construct form components
        radiantIntensityRay
        (
            const opticalDOM& dom,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label  lambda,
            const absorptionEmissionModel& absEmmModel_
        );


    // Destructor
    ~radiantIntensityRay();


    // Member functions

        // Edit

            //- Update optical intensity on i direction
            scalar correct();

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

            //- Add optical intensities from all the bands
            void addIntensity();


        // Access

            //- Return intensity
            inline const volScalarField& I() const;

            //- Return const access to the boundary heat flux
//            inline const volScalarField& Qr() const;

            //- Return non-const access to the boundary heat flux
//            inline volScalarField& Qr();

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the number of bands
            inline scalar nLambda() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            //- Return the optical intensity for a given wavelength
            inline const volScalarField& ILambda(const label lambdaI) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radiantIntensityRayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
