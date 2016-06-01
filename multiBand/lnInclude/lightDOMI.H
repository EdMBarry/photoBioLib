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

inline const Foam::light::lightIntensityRay&
Foam::light::lightDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}

/*
inline const Foam::volScalarField&
Foam::light::lightDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}
*/

inline Foam::label Foam::light::lightDOM::iRegion() const
{
    return iRegion_;
}

inline Foam::label Foam::light::lightDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::light::lightDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::light::lightDOM::nAngle() const
{
    return nAngle_;
}


inline Foam::label Foam::light::lightDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::light::lightDOM::nBand() const
{
    return nBand_;
}


inline Foam::label Foam::light::lightDOM::NumPixelPhi() const
{
    return NumPixelPhi_;
}

inline Foam::label Foam::light::lightDOM::NumPixelTheta() const
{
    return NumPixelTheta_;
}


/*
inline const Foam::volScalarField& Foam::light::lightDOM::k() const            // changed  from a_
{
    return k_;
}


inline const Foam::volScalarField& Foam::light::lightDOM::kLambda              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return kLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::light::lightDOM::sLambda              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return sLambda_[lambdaI];
}

*/

inline  Foam::scalar Foam::light::lightDOM::a              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return aBand_[lambdaI];
}


inline  Foam::scalar Foam::light::lightDOM::s              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return sBand_[lambdaI];
}

inline  Foam::scalar Foam::light::lightDOM::pf0              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return pf0_[lambdaI];
}

inline const Foam::volScalarField& Foam::light::lightDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::light::lightDOM::diffusionScatter() const
{
    return diffusionScatter_;
}




/*
inline const Foam::volScalarField& Foam::light::lightDOM::Qr() const
{
    return Qr_;
}

*/



// ************************************************************************* //
