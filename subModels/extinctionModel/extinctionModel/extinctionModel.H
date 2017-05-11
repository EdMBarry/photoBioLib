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
    Model to supply absorption and emission coefficients for photoBio
    modelling

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
 //       const fvMesh& mesh_;


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
                const dictionary& dict
  // ,            const fvMesh& mesh
            ),
            (dict)
   //         (dict, mesh)
        );


    // Constructors

        //- Construct from components
        extinctionModel
        (
            const dictionary& dict
  //          , const fvMesh& mesh
        );


    //- Selector
    static autoPtr<extinctionModel> New
    (
        const dictionary& dict
      //  ,const fvMesh& mesh
    );


    //- Destructor
    virtual ~extinctionModel();


    // Member Functions

        // Access

            //- Reference to the mesh
  /*          inline const fvMesh& mesh() const
            {
                return mesh_;
            }
*/
            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }


            // Extinction coefficient

   //         virtual tmp<volScalarField> k( label bandI ) const;

            virtual scalar k(const label bandI ) const;
 
            virtual scalar a(const label bandI ) const;
                       
            virtual scalar s(const label bandI ) const;
            
            // Absorption coefficient
  //          virtual tmp<volScalarField> a(const label bandI = 0) const;
            
            // Scattering coefficient
 //           virtual tmp<volScalarField> s( label bandI ) const;

 /*           //- Const access to the number of bands - defaults to 1 for grey
            //  Extinction
//            virtual label nBands() const;

            //- Const access to the bands - defaults to Vector2D::one for grey
            //  Extinction
 //           virtual const Vector2D<scalar>& bands(const label n) const;

            //- Flag for whether the Extinction/emission is for a grey gas
 //           virtual bool isGrey() const;

            //- Add radiative intensity for ray i
            virtual tmp<volScalarField> addIntensity
            (
                const label rayI,
                const volScalarField& ILambda
            ) const;

            //- Correct Extinction coefficients
            virtual void correct
            (
                volScalarField& a,
                PtrList<volScalarField>& aj
            ) const;
            
*/

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

