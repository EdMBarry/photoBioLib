/*---------------------------------------------------------------------------*\
  =========   	              |
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

#include "photoBioDOM.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "phaseFunctionModel.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(photoBioDOM, 0);
        addToRunTimeSelectionTable
        (
            photoBioModel,
            photoBioDOM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::photoBioDOM::photoBioDOM(const volScalarField& Irr)
:
    photoBioModel(typeName, Irr),
    iDomain_(readLabel(coeffs_.lookup("iDomain"))),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    diffusionScatter_
    (
         IOobject
        (
            "diffusionScatter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("diffusionScatter",dimMass/pow3(dimTime), 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nAngle_(0),
    nRay_(0),
    nBand_(coeffs_.lookupOrDefault<label>("nBand", 1)),
    NumPixelPhi_(coeffs_.lookupOrDefault<label>("NumPixelPhi", 1)),
    NumPixelTheta_(coeffs_.lookupOrDefault<label>("NumPixelTheta", 1)),
    GLambda_(nBand_),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50))
{
    Info<< "photoBioDOM number of Bands " << nBand_ << endl;

    // 3D
    if (mesh_.nSolutionD() == 3)
    {
		nAngle_ = 4*nPhi_*nTheta_;
        nRay_ = nAngle_*nBand_;
        IRay_.setSize(nRay_);
        scalar deltaPhi   =  pi /(2.0*nPhi_);
        scalar deltaTheta =  pi  /nTheta_;
        label i = 0;
        for (label iBand = 0; iBand < nBand_; iBand++)
        {
            for (label n = 1; n <= nTheta_; n++)
            {
                for (label m = 1; m <= 4*nPhi_; m++)
                {
    				label iAngle = m-1 + (n-1)*4*nPhi_;
                    scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                    scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                    IRay_.set
                    (
                        i,
                        new photoBioIntensityRay
                        (
                            *this,
                            mesh_,
                            iBand,
                            iAngle,
                            phii,
                            thetai,
                            deltaPhi,
                            deltaTheta
                        )
                    );
                    i++;
                }
            }
		}
    }
    else if (mesh_.nSolutionD() == 2)    //2D (X & Y)
    {
        if (mesh_.solutionD()[vector::Z] != -1)
        {
            FatalErrorInFunction
                << "Currently 2D solution is limited to the x-y plane"
                << exit(FatalError);
        }
        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        nAngle_ = 4*nPhi_;
        nRay_ = nAngle_*nBand_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi /(2.0*nPhi_);
        label i = 0;
        for (label iBand = 0; iBand < nBand_; iBand++)
		{
			for (label iAngle = 0; iAngle < 4*nPhi_; iAngle++)
			{
            scalar phii = (2.0*iAngle + 1.0)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new photoBioIntensityRay
                (
                    *this,
                    mesh_,
                    iBand,
                    iAngle,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta
                )
            );
            i++;
			}
		}
    }
    else    //1D (X)
    {
        if (mesh_.solutionD()[vector::X] != 1)
        {
            FatalErrorInFunction
                << "Currently 1D solution is limited to the x-direction"
                << exit(FatalError);
        }
        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        nAngle_ = 2;
        nRay_ = nAngle_*nBand_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi;
        label i = 0;
        for (label iBand = 0; iBand < nBand_; iBand++)
        {
            for (label m = 1; m <= 2; m++)
            {
			    label iAngle = m-1;
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new photoBioIntensityRay
                    (
                        *this,
                        mesh_,
                        iBand,
                        iAngle,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta
                    )
                );
                i++;
            }
        }
    }

    Info<< "photoBioDOM : Allocated " << IRay_.size() << nl;

   	sBand_.setSize(nBand_);
   	aBand_.setSize(nBand_);
    forAll(aBand_, iBand)
    {
		aBand_[iBand] = extinction_->a(iBand);
		sBand_[iBand] = extinction_->s(iBand);
	}


    forAll(GLambda_, iBand)
    {
        GLambda_.set
        (
            iBand,
            new volScalarField
            (
                IOobject
                (
                    "GLambda_" + Foam::name(iBand) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                G_                      // changed  from a_
            )
        );
    }

    phaseFunctionModel_ =  phaseFunctionModel::New(*this,coeffs_, mesh_.nSolutionD());

	if(phaseFunctionModel_->inScatter())
    {
        pf0_.setSize(nBand_ * nAngle_);
        for(label iBand = 0 ; iBand < nBand_; iBand++)
	    {
            for( label iAngle = 0 ; iAngle < nAngle_; iAngle++)
            {
                pf0_[iAngle+iBand*nAngle_] = phaseFunctionModel_->correct(iAngle,iAngle,iBand);
            }
        }
	}

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::photoBio::photoBioDOM::~photoBioDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::photoBio::photoBioDOM::read()
{
    if (photoBioModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::photoBio::photoBioDOM::calculate()
{
    scalar maxResidual = 0.0;
    label rayJ;
	label iBand = 0;
    label iAngle = 0;
    label radIter = 0;
    scalar maxBandResidual = 0.0;

    do
    {
        radIter++;
        forAll(IRay_, rayI)  //		label	rayI = 0;
		{
			maxResidual = 0.0;
			iBand = IRay_[rayI].iBand();
			iAngle = IRay_[rayI].iAngle();

            Info << "photoBio solver iDomain:  " << iDomain_ << "    iter: " << radIter << "    iBand: "<< iBand << "    iAngle  : "<< iAngle << endl;
            Info<< endl;

			if(phaseFunctionModel_->inScatter())
            {
                diffusionScatter_ = dimensionedScalar("diffusionScatter",dimMass/pow3(dimTime), 0.0) ;

			    for (label jAngle = 0; jAngle < nAngle_; jAngle++)
				{
					rayJ = jAngle + iBand*nAngle_;
					if(rayI != rayJ )       //rayCos = 1; 	//if(rayCos < -1) rayCos = -1;
                    {
			     	    scalar rayCos = IRay_[rayI].d() & IRay_[rayJ].d();
                        if(rayCos > 0 )
                        {
		                    diffusionScatter_ = diffusionScatter_ + IRay_[rayJ].I()*phaseFunctionModel_->correct(rayI,rayJ, iBand)*IRay_[rayJ].omega();
                        }
                    }
				}
			}

            IRay_[rayI].updateBoundary();
            maxBandResidual = IRay_[rayI].correct();            // solve the RTE equation
			maxResidual = max(maxBandResidual, maxResidual);
        }
    } while(maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}


void Foam::photoBio::photoBioDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    label rayI;
    forAll(GLambda_, iBand)
    {
	//		Info << "iBand: " << iBand << endl;
        GLambda_[iBand] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
        for (label iAngle = 0; iAngle < nAngle_; iAngle++)
	    {
            rayI = iAngle + iBand*nAngle_;
            GLambda_[iBand] += IRay_[rayI].I()*IRay_[rayI].omega();  // convert the radiance to irradiance
        }
			G_ += GLambda_[iBand];
	}
}


void Foam::photoBio::photoBioDOM::setRayId
(
    const word& name,
    label& rayId
) const
{
    // assuming name is in the form: CHARS_iBand_iAngle

    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    label ib = readLabel(IStringStream(name.substr(i1+1, i2-1))());

    label ia = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());

    rayId = nAngle_*ib + ia;

}

void Foam::photoBio::photoBioDOM::dirToRayId
(
    const vector& dir,
    const  label& iBand,
    label& rayId
) const
{
	scalar tTheta = Foam::acos(dir.z()/mag(dir));
	scalar tPhi;
	if(dir.x() != 0 )
	{
	tPhi = Foam::atan(dir.y()/dir.x());
	if(dir.x() < 0 && dir.y() > 0 ) tPhi = tPhi + pi;
	if(dir.x() < 0 && dir.y() < 0 ) tPhi = tPhi + pi;
	if(dir.x() > 0 && dir.y() < 0 ) tPhi = tPhi + 2*pi;
	}
	else
	{
		if(dir.y() > 0 ) tPhi = pi/2.0;
		if(dir.y() < 0 ) tPhi = 3*pi/2.0;;
	}


    label iPhi = label(tPhi/deltaPhi);
    label iTheta = label(tTheta/deltaTheta);

    rayId = nAngle_*iBand + iTheta*4*nPhi_ + iPhi;



}



/*
    // Construct extinction field for each wavelength
    forAll(kLambda_, iBand)     // changed  from aLambda
    {
        kLambda_.set
        (
            iBand,
            new volScalarField
            (
                IOobject
                (
                    "kLambda_" + Foam::name(iBand) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                k_                      // changed  from a_
            )
        );
    }

    forAll(sLambda_, iBand)     // changed  from aLambda
    {
        sLambda_.set
        (
            iBand,
            new volScalarField
            (
                IOobject
                (
                    "sLambda_" + Foam::name(iBand) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                k_                      // changed  from a_
            )
        );
    }
*/
// ************************************************************************* //
