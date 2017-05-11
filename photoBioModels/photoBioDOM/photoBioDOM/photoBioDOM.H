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
    Foam::photoBio::photoBioDOM

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter.

    Available absorptionEmission models:
        greyMeanAbsoprtionEmission
        wideBandabsorptionEmissionEmission

    i.e. dictionary
    photoBioDOMCoeffs
    {
        nPhi    1;          // azimuthal angles in PI/2 on X-Y.(from Y to X)
        nTheta  2;          // polar angles in PI (from Z to X-Y plane)
        convergence 1e-4;   // convergence criteria for photoBio iteration
    }

    solverFreq   1; // Number of flow iterations per photoBio iteration

    The total number of solid angles is  4*nPhi*nTheta.

    In 1D the direction of the rays is X (nPhi and nTheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)

SourceFiles
    photoBioDOM.C

\*---------------------------------------------------------------------------*/

#ifndef photoBioModelphotoBioDOM_H
#define photoBioModelphotoBioDOM_H

#include "photoBioModel.H"
#include "photoBioIntensityRay.H"
#include "phaseFunctionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

/*---------------------------------------------------------------------------*\
                           Class photoBioDOM Declaration
\*---------------------------------------------------------------------------*/

class photoBioDOM
:
    public photoBioModel
{
    // Private data

        //- Incident photoBio  [lum]
        label iRegion_ ;

        volScalarField G_;

        volScalarField diffusionScatter_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        // - Number of solid angles
        label nAngle_;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nBand_;

        label iDomain_;

        label NumPixelPhi_;

        label NumPixelTheta_;

        PtrList<volScalarField> GLambda_;

        DynamicList<scalar, 1, 0>   aBand_;

        DynamicList<scalar, 1, 0>   sBand_;

        DynamicList<scalar, 1, 0>   pf0_;

           //- List of pointers to photoBio intensity rays
        PtrList<photoBioIntensityRay> IRay_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        scalar deltaPhi   ;

        scalar deltaTheta ;

        autoPtr<phaseFunctionModel> phaseFunctionModel_;



    // Private member functions

        //- Disallow default bitwise copy construct
        photoBioDOM(const photoBioDOM&);

        //- Disallow default bitwise assignment
        void operator=(const photoBioDOM&);



public:

    //- Runtime type information
    TypeName("photoBioDOM");


    // Constructors

        //- Construct from components
        photoBioDOM(const volScalarField& intensity);


    //- Destructor
    virtual ~photoBioDOM();


    // Member functions

        // Edit

            //- Solve photoBio equation(s)
            void calculate();

            //- Read photoBio properties dictionary
            bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayId
            (
                const word& name,
                label& rayId
            ) const;


            void dirToRayId
            (
                const vector& dir,
                const  label& iBand,
                label& rayId
            ) const;


            //- Ray intensity for rayI
            inline const photoBioIntensityRay& IRay(const label rayI) const;

            inline label iRegion() const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of angles in phi
            inline label nAngle() const;


            //- Number of rays
            inline label nRay() const;

            //- Number of wavelengths
            inline label nBand() const;

            inline label NumPixelPhi() const;

            inline label NumPixelTheta() const;

            inline  scalar a(const label lambdaI) const;

            inline  scalar s(const label lambdaI) const;

            inline  scalar pf0(const label lambdaI) const;

            //- Const access to incident photoBio field
            inline const volScalarField& G() const;

            inline const volScalarField& diffusionScatter() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "photoBioDOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
