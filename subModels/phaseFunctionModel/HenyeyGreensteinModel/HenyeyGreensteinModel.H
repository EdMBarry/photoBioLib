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
    Foam::photoBio::constantScatter

Description
    Constant photoBio scatter coefficient

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef HenyeyGreensteinModel_H
#define HenyeyGreensteinModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "scalar.H"

#include "phaseFunctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

class photoBioDOM;
/*---------------------------------------------------------------------------*\
                           Class constantScatter Declaration
\*---------------------------------------------------------------------------*/

class HenyeyGreensteinModel : public phaseFunctionModel
{
private:

    // Private data
    Switch inScatter_;
    dictionary coeffsDict_;
    label subAngleNum_ ;
    DynamicList<scalar, 1, 0> g_;
    DynamicList<scalar, 1, 0> phaseFunction_;
    
public:

    //- Runtime type information
    TypeName("HenyeyGreensteinModel");

    // Constructors
    //- Construct from components
    HenyeyGreensteinModel
    (
        const photoBioDOM& dom,
        const dictionary& dict,
        const label& nDim
    );


    // Destructor
    virtual ~HenyeyGreensteinModel();


    // Member Operators
    Foam::scalar correct
    (
        const label rayI,
        const label rayJ,
        const label iBand
    ) const;
        
    Foam::scalar hg3d
    (
        const scalar cosV,
        const scalar g
    ) const;
    
    Foam::scalar hg2d
    (
        const scalar cosV,
        const scalar g
    ) const;
    
};
    
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
