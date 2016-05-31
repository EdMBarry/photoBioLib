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
    Foam::optical::P1

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

SourceFiles
    P1.C

\*---------------------------------------------------------------------------*/

#ifndef opticalModelP1_H
#define opticalModelP1_H

#include "opticalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                           Class P1 Declaration
\*---------------------------------------------------------------------------*/

class P1
:
    public opticalModel
{
    // Private data

        //- Incident optical / [W/m2]
        volScalarField G_;

        //- Total radiative heat flux [W/m2]
   //     volScalarField Qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
   //     volScalarField e_;

        //- Emission contribution
   //     volScalarField E_;


    // Private member functions

        //- Disallow default bitwise copy construct
        P1(const P1&);

        //- Disallow default bitwise assignment
        void operator=(const P1&);


public:

    //- Runtime type information
    TypeName("opticalP1");


    // Constructors

        //- Construct from components
        P1(const volScalarField& intensity);


    // Destructor
    virtual ~P1();


    // Member functions

        // Edit

            //- Solve optical equation(s)
            void calculate();

            //- Read optical properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
   //         virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
   //         virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
