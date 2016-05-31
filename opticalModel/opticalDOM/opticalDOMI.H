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

\*---------------------------------------------------------------------------*/

inline const Foam::optical::radiantIntensityRay&
Foam::optical::opticalDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::optical::opticalDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::optical::opticalDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::optical::opticalDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::optical::opticalDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::optical::opticalDOM::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::optical::opticalDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::optical::opticalDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::optical::opticalDOM::G() const
{
    return G_;
}

/*
inline const Foam::volScalarField& Foam::optical::opticalDOM::Qr() const
{
    return Qr_;
}

*/



// ************************************************************************* //
