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
    Foam::optical::opticalDOM

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter.

    Available absorptionEmission models:
        greyMeanAbsoprtionEmission
        wideBandabsorptionEmissionEmission

    i.e. dictionary
    opticalDOMCoeffs
    {
        nPhi    1;          // azimuthal angles in PI/2 on X-Y.(from Y to X)
        nTheta  2;          // polar angles in PI (from Z to X-Y plane)
        convergence 1e-4;   // convergence criteria for optical iteration
    }

    solverFreq   1; // Number of flow iterations per optical iteration

    The total number of solid angles is  4*nPhi*nTheta.

    In 1D the direction of the rays is X (nPhi and nTheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)

SourceFiles
    opticalDOM.C

\*---------------------------------------------------------------------------*/

#ifndef opticalModelopticalDOM_H
#define opticalModelopticalDOM_H

#include "opticalModel.H"
#include "opticalIntensityRay.H"
#include "phaseFunctionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                           Class opticalDOM Declaration
\*---------------------------------------------------------------------------*/

class opticalDOM
:
    public opticalModel
{
    // Private data

        //- Incident optical  [lum]
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

         label NumPixelPhi_;
         
         label NumPixelTheta_;


        PtrList<volScalarField> GLambda_; 
         
        DynamicList<scalar, 1, 0>   aBand_;
        
        DynamicList<scalar, 1, 0>   sBand_;
        
        DynamicList<scalar, 1, 0>   pf0_;
         
           //- List of pointers to optical intensity rays
        PtrList<opticalIntensityRay> IRay_;
   
        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;
       
        scalar deltaPhi   ;
        
        scalar deltaTheta ;

        autoPtr<phaseFunctionModel> phaseFunctionModel_;
         

        
    // Private member functions

        //- Disallow default bitwise copy construct
        opticalDOM(const opticalDOM&);

        //- Disallow default bitwise assignment
        void operator=(const opticalDOM&);



public:

    //- Runtime type information
    TypeName("opticalDOM");


    // Constructors

        //- Construct from components
        opticalDOM(const volScalarField& intensity);


    //- Destructor
    virtual ~opticalDOM();


    // Member functions

        // Edit

            //- Solve optical equation(s)
            void calculate();

            //- Read optical properties dictionary
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
            inline const opticalIntensityRay& IRay(const label rayI) const;

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
             
            //- Const access to incident optical field
            inline const volScalarField& G() const;
            
            inline const volScalarField& diffusionScatter() const;
          

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "opticalDOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
