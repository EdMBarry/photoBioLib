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

#include "lightIntensityRay.H"
#include "fvm.H"
#include "lightDOM.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::light::lightIntensityRay::rayId(0);

const Foam::word
Foam::light::lightIntensityRay::intensityPrefix("I");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::lightIntensityRay::lightIntensityRay
(
    const lightDOM& dom,
    const fvMesh&   mesh,
    const label     iBand,
    const label     iAngle,
    const scalar    phi,
    const scalar    theta,
    const scalar    deltaPhi,
    const scalar    deltaTheta
)
:
    dom_(dom),
    mesh_(mesh),
    d_(vector::zero),
    dAve_(vector::zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    iBand_(iBand),
    iAngle_(iAngle)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

     
   omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    //  omega_ = deltaPhi;
     
    d_ = vector(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
   //   d_ = vector(cosPhi, sinPhi,0.0);  
      
   dAve_ = vector
    (
       cosPhi 
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
       
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
       
        0.5*deltaPhi
        *Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );

 /*
    dAve_ = vector(2*cosPhi*Foam::sin(0.5*deltaPhi),2*sinPhi*Foam::sin(0.5*deltaPhi),0.0);
  */
      autoPtr<volScalarField> IDefaultPtr;


        IOobject IHeader
        (
            intensityPrefix + "_" + name(iBand)+ "_" + name(iAngle),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        // check if field exists and can be read
        if (IHeader.headerOk())
        {
           I_.set(new volScalarField(IHeader, mesh_));
        }
        else
        {
			
            // Demand driven load the IDefault field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                           "IDefault",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }

            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

			I_.set(new volScalarField(noReadHeader, IDefaultPtr()));
			
        }
    
    
    rayId++;
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::lightIntensityRay::~lightIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::scalar Foam::light::lightIntensityRay::correct()  //(const volScalarField& inScTerm)
{
    // reset boundary heat flux to zero
    // Qr_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
    //    Qr_.boundaryField() = 0.0;

    scalar maxResidual = -GREAT;
    scalar eqnResidual ;
      
      scalar selfScat = 0.0;
      if( dom_.s(iBand_) != 0.0 )  selfScat =  dom_.pf0(iAngle_ + iBand_*dom_.nAngle()) ;
     
     dimensionedScalar k("k", dimless/dimLength, dom_.a(iBand_)+dom_.s(iBand_)*(1.0- selfScat*omega_) );  // 0.5
     dimensionedScalar s("s", dimless/dimLength, dom_.s(iBand_)); 
  
     
     const volScalarField& ds = dom_.diffusionScatter();
     const label npP = dom_.NumPixelPhi();
     label npT  = 1;
     if(mesh_.nSolutionD() == 3) npT = dom_.NumPixelTheta();
     
      vector tV = vector(0,0,0);                 
	  surfaceScalarField Ji(tV & mesh_.Sf());
	  
	   const scalar  dp = pi /(2.0*dom_.nPhi())/npP ; 
	   const scalar  dt = pi /dom_.nTheta()/npT;

      for(label m=0; m < npT;m++)
      {
		  for(label n=0; n < npP;n++)
		  {
			  	scalar nP = (2.0*m -npP +1.0)*dp/2.0 + phi_; 
				scalar nT = (2.0*n -npT +1.0)*dt/2.0 + theta_;
				vector nD = vector (   
				Foam::cos(nP)*Foam::sin(0.5*dp)*(dt - Foam::cos(2.0*nT)*Foam::sin(dt)),
				Foam::sin(nP)*Foam::sin(0.5*dp)*(dt - Foam::cos(2.0*nT)*Foam::sin(dt)),
				0.5*dp*Foam::sin(2.0*nT)*Foam::sin(dt));
				Ji = Ji + (nD & mesh_.Sf());
			}
		}
				

       fvScalarMatrix IiEq
        (
            fvm::div(Ji, I_(), "div(Ji,Ii_h)")
          + fvm::Sp(k*omega_, I_())            
          == s*ds*omega_  
        );

        IiEq.relax();

     eqnResidual = solve(IiEq,mesh_.solver("Ii")).initialResidual();

     maxResidual = max(eqnResidual, maxResidual);
 
     return maxResidual;
    
    
}



void Foam::light::lightIntensityRay::updateBoundary()  
{

           I_->correctBoundaryConditions();

}
// ************************************************************************* //
