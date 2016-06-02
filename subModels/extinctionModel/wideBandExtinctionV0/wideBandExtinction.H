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
    Foam::optical::bioMediaAbsorption

Description
    Constant optical absorption and emission coefficients for continuous
    phase

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef opticalWideBandExtinction_H
#define opticalWideBandExtinction_H

// #include "HashTable.H"
#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                  Class wideBandExtinction Declaration
\*---------------------------------------------------------------------------*/

class wideBandExtinction
:
    public extinctionModel
{
	
    // Private data
    
        //- Maximum number of bands
        
        //- Extinction model dictionary
        dictionary coeffsDict_;

       //- Extinction coefficient / [1/m]       //     dimensionedScalar a_;
       
        DynamicList<scalar, 1, 0>   kBand_;
               
        DynamicList<scalar, 1, 0>   aBand_;
        
        DynamicList<scalar, 1, 0>   sBand_;
        
        //- Bands
        label nBands_;
        
        //- Bands
        DynamicList<Vector2D<scalar>, 1,0> iBands_;

        //- Proportion of the heat released rate emitted
   //     FixedList<scalar, maxBands_> iEhrrCoeffs_;
        
        // Total wave length covered by the bands
        scalar totalWaveLength_;
        
  //      mutable interpolationLookUpTable<scalar> lookUpTable_;

        
public:

    //- Runtime type information
    TypeName("wideBandExtinction");


    // Constructors

        //- Construct from components
        wideBandExtinction
        (
            const dictionary& dict
//            , const fvMesh& mesh
        );


    // Destructor
    virtual ~wideBandExtinction();


    // Member Operators

        // Access

        //- Extinction coefficient for continuous phase
          scalar k( const label bandI ) const;
                 
          scalar a( const label bandI ) const;

		  scalar s( const label bandI ) const;

 //       tmp<volScalarField> s( label bandI ) const;
         
         
        inline const Vector2D<scalar>& bands(const label i) const
        {
              return iBands_[i];
        }         
 /*        
        inline bool isGrey() const
        {
              return false;                // this needs to be in the optical properties file
        }

        //- Number of bands
        inline label nBands() const
        {
              return nBands_;
        }
*/
         //- Lower and upper limit of band i


        //- Add contribution of ILambda to the total radiative intensity in
        //  direction i
 /*       tmp<volScalarField> addIntensity
        (
              const label i,
              const volScalarField& ILambda
        ) const;


        void correct
        (
              volScalarField& k_,
              PtrList<volScalarField>& kLambda
        ) const;

*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
