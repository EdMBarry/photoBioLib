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
    Foam::light::constantScatter

Description
    Constant light scatter coefficient

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef userDefinedPhaseFunction_H
#define userDefinedPhaseFunction_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "scalar.H"

#include "phaseFunctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

class lightDOM;

/*---------------------------------------------------------------------------*\
                           Class constantScatter Declaration
\*---------------------------------------------------------------------------*/

class userDefinedPhaseFunction
:
    public phaseFunctionModel
{
public:


	    DynamicList<scalar, 1, 0>  phaseFunction_;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;
               
        Switch inScatter_;
        
        label userPhaseFuncAngleNum_ ;
        label subAngleNum_ ;
        DynamicList<scalar, 1, 0>  userPhaseFunc_;

public:

    //- Runtime type information
    TypeName("userDefinedPhaseFunction");


    // Constructors

        //- Construct from components
        userDefinedPhaseFunction
        (
                const lightDOM& dom,
				const 	dictionary& dict,
				const label& nDim
        );


    // Destructor

       virtual ~userDefinedPhaseFunction();


    // Member Operators

       inline bool inScatter() const
        {
            return inScatter_;
        }
       

         Foam::scalar correct
        (
          const  label  rayI,
          const  label  rayJ,
            const label  iBand
        ) const;
        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
