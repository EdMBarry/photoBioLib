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

\*---------------------------------------------------------------------------*/

#include "userDefinedPhaseFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "lightDOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
using namespace Foam::constant::mathematical;
namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(userDefinedPhaseFunction, 0);

        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            userDefinedPhaseFunction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::userDefinedPhaseFunction::userDefinedPhaseFunction
(
                const lightDOM& dom,
				const dictionary& dict,
				const label& nDim
)
:
    phaseFunctionModel(dom,dict,nDim),
  //  coeffsDict_(dict.subDict(typeName + "Coeffs"))
    coeffsDict_(dict)
{
    const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    
    functionDicts.lookup("inScatter") >> inScatter_;
    
  	if(inScatter_)
    {  
       nBand_ = dom_.nBand();   
	   nAngle_ = dom_.nAngle();
	     
		     functionDicts.lookup("subAngleNum") >> subAngleNum_;
		     
    functionDicts.lookup("userPhaseFuncAngleNum") >> userPhaseFuncAngleNum_;
    
    userPhaseFunc_.setSize(nBand_*userPhaseFuncAngleNum_);   
     
	functionDicts.lookup("userPhaseFunction") >> userPhaseFunc_;
	   
     scalar deltaAngle = Foam::constant::mathematical::pi/scalar(userPhaseFuncAngleNum_);   
       
     const  label nPhi = dom_.nPhi();  
     const  label nTheta = dom_.nTheta();   
     const  scalar deltaPhi   =  pi /(2.0*nPhi);
     const  scalar deltaTheta =  pi  /nTheta;
     
     phaseFunction_.setSize(nAngle_*nAngle_*nBand_);    
     forAll(phaseFunction_,i)     phaseFunction_[i] = 0.0;
     
  //   	Info<< "phaseFunction_.setSize " << nAngle_*nAngle_*nBand_ << endl;
       
    if (nDim == 3)    //3D 
    {
	scalar  dp = deltaPhi/subAngleNum_;
	scalar  dt = deltaTheta/subAngleNum_;
	
	for(label iband = 0; iband < nBand_ ; iband++)
	{
		for(label i = 0; i<nAngle_ ; i++)
		{
			label rayI = i + iband*nAngle_;  
			scalar pfSum = 0;
			for(label j = 0; j<nAngle_ ; j++)
			{
				label  rayJ = j + iband*nAngle_;  
				label  idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
				for(label m = 0; m < subAngleNum_ ; m++)
				{
				for(label n = 0; n < subAngleNum_ ; n++)
				{
					scalar nP = (2.0*m -subAngleNum_ +1.0)*dp/2.0 + dom.IRay(rayJ).phi(); 
					scalar nT = (2.0*n -subAngleNum_ +1.0)*dt/2.0 + dom.IRay(rayJ).theta();
					scalar nOmega = 2*sin(nT)*sin(dt/2)*dp;
					vector nD = vector (sin(nT)*cos(nP), sin(nT)*sin(nP), cos(nT));
					scalar cosV = dom.IRay(rayI).d()  & nD;
					label  i0 =0;
					if(cosV > 1) i0=0;
					else if(cosV < -1) i0 = userPhaseFuncAngleNum_-1;
					else  i0 =label(Foam::acos(cosV)/deltaAngle+0.5);
					
					phaseFunction_[idx] = phaseFunction_[idx] + userPhaseFunc_[i0]*nOmega;
				}
				}
				pfSum = pfSum + phaseFunction_[idx];
				phaseFunction_[idx] = phaseFunction_[idx]/dom.IRay(rayI).omega() ;
			}
			
			for(label j = 0; j<nAngle_ ; j++)
			{
			label  idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
			phaseFunction_[idx] = phaseFunction_[idx]/pfSum;
			}
		}
	}
	}
     
    if (nDim == 2)    //2D 
    {
	scalar  dp = deltaPhi/subAngleNum_;
	scalar  dt = pi/90;
	
	for(label iband = 0; iband < nBand_ ; iband++)
	{
		for(label i = 0; i<nAngle_ ; i++)
		{
			label rayI = i + iband*nAngle_;  
			scalar pfSum = 0;
			for(label j = 0; j<nAngle_; j++)
			{
				label  rayJ = j + iband*nAngle_;  
				label  idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
				for(label m = 0; m < subAngleNum_ ; m++)
				{
				for(label n = 0; n < 90 ; n++)
				{
					scalar nP = (2.0*m -subAngleNum_ +1.0)*dp/2.0 + dom.IRay(rayJ).phi(); 
					scalar nT = (2.0*n + 1.0)*dt/2.0;
					scalar nOmega = 2*sin(nT)*sin(dt/2)*dp;
					vector nD = vector (sin(nT)*cos(nP), sin(nT)*sin(nP), cos(nT));
					scalar cosV = dom.IRay(rayI).d()  & nD;
					label  i0 = label(Foam::acos(cosV)/deltaAngle+0.5);
					phaseFunction_[idx] = phaseFunction_[idx] + userPhaseFunc_[i0]*sin(nT)*nOmega;
				}
				}
				
				pfSum = pfSum + phaseFunction_[idx];
				phaseFunction_[idx] = phaseFunction_[idx]/dom.IRay(rayI).omega() ;
			}
			
			for(scalar j = 0; j<nAngle_ ; j++)
			{
			label  idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
			phaseFunction_[idx] = phaseFunction_[idx]/pfSum;  // 4/pi
			}
		}
	}
	}
	
	
//	 for(label i = 0; i< 40 ; i++)  Info<< phaseFunction_[i]  << endl; 
	
	
	}
	
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::userDefinedPhaseFunction::~userDefinedPhaseFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




 Foam::scalar  Foam::light::userDefinedPhaseFunction::correct
(
          const   label rayI,
          const   label rayJ,
          const   label  iBand
) const
{
	  return phaseFunction_[rayJ + rayI*nAngle_ +iBand*nAngle_*nAngle_ ];
}


/*
 *  Foam::scalar  Foam::light::userDefinedPhaseFunction::correct
(
          const   label rayI,
          const   label rayJ,
          const   label  iBand
) const
{
	   scalar deltaAngle = Foam::constant::mathematical::pi/scalar(phaseFunctionAngleNum_-1);   
       scalar angle = 	Foam::acos(rayCos);
       label  i0 = label(angle/deltaAngle+0.2);
       
       if(i0 >= 0 && i0 < phaseFunctionAngleNum_)
       { 
               return  phaseFunc_[i0 + iBand*phaseFunctionAngleNum_];
			//	return 0.0;
		}
		else
		{ 
			return 0.0;
		}
       
}*/

// ************************************************************************* //
