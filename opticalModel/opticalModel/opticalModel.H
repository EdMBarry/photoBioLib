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

Namespace
    Foam::optical

Description
    Namespace for optical modelling

Class
    Foam::optical::opticalModel

Description
    Top level model for optical modelling

SourceFiles
    opticalModel.C

\*---------------------------------------------------------------------------*/

#ifndef opticalModel_H
#define opticalModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

// Forward declaration of classes
class absorptionEmissionModel;
class scatterModel;

/*---------------------------------------------------------------------------*\
                       Class opticalModel Declaration
\*---------------------------------------------------------------------------*/

class opticalModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;


        //- Model specific dictionary input parameters
        Switch optical_;

        //- optical model dictionary
        dictionary coeffs_;

        //- optical solver frequency - number flow solver iterations per
        //  optical solver iteration
        label solverFreq_;


        // References to the optical sub-models

            //- Absorption/emission model
            autoPtr<absorptionEmissionModel> absorptionEmission_;

            //- Scatter model
            autoPtr<scatterModel> scatter_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        opticalModel(const opticalModel&);

        //- Disallow default bitwise assignment
        void operator=(const opticalModel&);


public:

    //- Runtime type information
    TypeName("opticalModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             opticalModel,
             dictionary,
             (
                 const volScalarField& intensity
             ),
             (intensity)
         );


    // Constructors

        //- Null constructor
        opticalModel(const volScalarField& intensity);

        //- Construct from components
        opticalModel(const word& type, const volScalarField& intensity);


    // Selectors

         //- Return a reference to the selected optical model
         static autoPtr<opticalModel> New(const volScalarField& intensity);


    //- Destructor
    virtual ~opticalModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve optical equation(s)
            virtual void calculate() = 0;

            //- Read opticalProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Source term component (for power of T^4)
  //          virtual tmp<volScalarField> Rp() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
