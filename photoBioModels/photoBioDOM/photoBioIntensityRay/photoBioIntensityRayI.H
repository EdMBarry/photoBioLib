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


 inline const Foam::volScalarField& Foam::photoBio::photoBioIntensityRay::I() const
{
    return I_;
}


inline const Foam::vector& Foam::photoBio::photoBioIntensityRay::d() const
{
    return d_;
}


inline const Foam::vector& Foam::photoBio::photoBioIntensityRay::dAve() const
{
    return dAve_;
}


inline Foam::label Foam::photoBio::photoBioIntensityRay::iBand() const
{
    return iBand_;
}


inline Foam::label Foam::photoBio::photoBioIntensityRay::iAngle() const
{
    return iAngle_;
}


inline Foam::scalar Foam::photoBio::photoBioIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::photoBio::photoBioIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::photoBio::photoBioIntensityRay::omega() const
{
    return omega_;
}

// ************************************************************************* //
