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

\*---------------------------------------------------------------------------*/

#include "sunOpticalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "opticalDOM.H"
#include "surfaceFields.H"

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sunOpticalFvPatchScalarField::
sunOpticalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    d0_(vector::zero),
    I0_(0)
{}


Foam::sunOpticalFvPatchScalarField::
sunOpticalFvPatchScalarField
(
    const sunOpticalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    d0_((ptf.d0_),
    I0_(ptf.I0_)
{}


Foam::sunOpticalFvPatchScalarField::
sunOpticalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    d0_(dict.lookup("origin")),
    I0_(readScalar(dict.lookup("omega")))
{
    // Evaluate the wall velocity
    updateCoeffs();
}


Foam::sunOpticalFvPatchScalarField::
sunOpticalFvPatchScalarField
(
    const sunOpticalFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchField<vector>(pivpvf),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_)
{}


Foam::sunOpticalFvPatchScalarField::
sunOpticalFvPatchScalarField
(
    const sunOpticalFvPatchScalarField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(pivpvf, iF),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sunOpticalFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate the rotating wall velocity from the specification of the motion
    const vectorField Up
    (
        (-omega_)*((patch().Cf() - origin_) ^ (axis_/mag(axis_)))
    );

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    const vectorField n(patch().nf());
    vectorField::operator=(Up - n*(n & Up));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::sunOpticalFvPatchScalarField::write(Ostream& os) const
{
    FvPatchScalarField::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << omega_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        FvPatchScalarField,
        sunOpticalFvPatchScalarField
    );
}

// ************************************************************************* //
