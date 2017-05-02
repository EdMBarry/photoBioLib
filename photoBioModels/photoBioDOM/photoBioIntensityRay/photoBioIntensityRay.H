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

        //-  photoBio ray intensity / [W/m2]

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

        //- Number of wavelengths/bands
        label iBand_;

        label iAngle_;

        //- List of pointers to photoBio intensity fields for given wavelengths
    //    PtrList<volScalarField> ILambda_;

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
            scalar correct();  //(const volScalarField& inScTerm)

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

            //- Add photoBio intensities from all the bands
  //          void addIntensity();


        // Access

            //- Return intensity
             inline  const volScalarField&  I() const;
         //        inline const volScalarField& I() const;

            //- Return const access to the boundary heat flux
//            inline const volScalarField& Qr() const;

            //- Return non-const access to the boundary heat flux
//            inline volScalarField& Qr();

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the total number of bands
    //        inline scalar nBand() const;

             //- Return the number of band
            inline  label iBand() const;

            inline label iAngle() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            const  Foam::surfaceScalarField JiCal() const;

            //- Return the photoBio intensity for a given wavelength
     //       inline const volScalarField& ILambda(const label lambdaI) const;

           void updateBoundary() ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "photoBioIntensityRayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
