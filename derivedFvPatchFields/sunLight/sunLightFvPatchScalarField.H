/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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
    Foam::sunOpticalFvPatchScalarField

Description
    Foam::sunOpticalFvPatchScalarField

SourceFiles
    sunOpticalFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sunOpticalFvPatchScalarField_H
#define sunOpticalFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sunOpticalFvPatch Declaration
\*---------------------------------------------------------------------------*/

class sunOpticalFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        vector  d0_ ;         // incident optical direction
        
        scalar  I0_;              // incident optical intensity


public:

    //- Runtime type information
    TypeName("sunOptical");


    // Constructors

        //- Construct from patch and internal field
        sunOpticalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sunOpticalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given sunOpticalFvPatchScalarField
        //  onto a new patch
        sunOpticalFvPatchScalarField
        (
            const sunOpticalFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sunOpticalFvPatchScalarField
        (
            const sunOpticalFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<FvPatchScalarField> clone() const
        {
            return tmp<FvPatchScalarField>
            (
                new sunOpticalFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sunOpticalFvPatchScalarField
        (
            const sunOpticalFvPatchScalarField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<FvPatchScalarField> clone
        (
             const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<FvPatchScalarField>
            (
                new sunOpticalFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        // Access functions


            //- Return the rotational speed
            scalar d0() const
            {
                return d0_;
            }

            //- Return non-const access to the origin of the rotation
            vector& I0()
            {
                return I0_;
            }



        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
