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
    Foam::photoBio::wideBandConstantExtinction

Description
    Model for constant absorption and scattering coefficients

SourceFiles
    wideBandConstantExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef photoBioWideBandExtinction_H
#define photoBioWideBandExtinction_H

#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

/*---------------------------------------------------------------------------*\
                  Class wideBandConstantExtinction Declaration
\*---------------------------------------------------------------------------*/

class wideBandConstantExtinction
:
    public extinctionModel
{
	
    // Private data
        
        //- Extinction model dictionary
        dictionary coeffsDict_;
        
        //- Switch for absorption model
        Switch absorption_;
        
        //- Switch of scattering model
        Switch scattering_;   
        
        //- Number of bands
        label nBands_;

        //- Absorption coefficients (one per band)
        DynamicList<scalar, 1, 0> ABand_;
        
        //- Scattering coefficients (one per band)
        DynamicList<scalar, 1, 0> SBand_;
       
public:

    //- Runtime type information
    TypeName("wideBandConstantExtinction");


    // Constructors

        //- Construct from components
        wideBandConstantExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    // Destructor
    virtual ~wideBandConstantExtinction();

    // Member Operators

        //- Correct the extinction coefficients
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
