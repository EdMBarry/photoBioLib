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

inline const Foam::photoBio::photoBioIntensityRay&
Foam::photoBio::photoBioDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


/*inline const Foam::volScalarField&		// previously commented out
Foam::photoBio::photoBioDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}
*/

inline Foam::label Foam::photoBio::photoBioDOM::iRegion() const
{
    return iRegion_;
}

inline Foam::label Foam::photoBio::photoBioDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::photoBio::photoBioDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::photoBio::photoBioDOM::nAngle() const
{
    return nAngle_;
}


inline Foam::label Foam::photoBio::photoBioDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::photoBio::photoBioDOM::nBand() const
{
    return nBand_;
}


inline Foam::label Foam::photoBio::photoBioDOM::NumPixelPhi() const
{
    return NumPixelPhi_;
}

inline Foam::label Foam::photoBio::photoBioDOM::NumPixelTheta() const
{
    return NumPixelTheta_;
}


/*
inline const Foam::volScalarField& Foam::photoBio::photoBioDOM::k() const            // changed  from a_
{
    return k_;
}


inline const Foam::volScalarField& Foam::photoBio::photoBioDOM::kLambda              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return kLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::photoBio::photoBioDOM::sLambda              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return sLambda_[lambdaI];
}

*/

inline  Foam::scalar Foam::photoBio::photoBioDOM::a              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return aBand_[lambdaI];
}


inline  Foam::scalar Foam::photoBio::photoBioDOM::s              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return sBand_[lambdaI];
}

inline  Foam::scalar Foam::photoBio::photoBioDOM::pf0              // changed  from kLambda_
(
    const label lambdaI
) const
{
    return pf0_[lambdaI];
}

inline const Foam::volScalarField& Foam::photoBio::photoBioDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::photoBio::photoBioDOM::diffusionScatter() const
{
    return diffusionScatter_;
}




/*
inline const Foam::volScalarField& Foam::photoBio::photoBioDOM::Qr() const
{
    return Qr_;
}

*/



// ************************************************************************* //
