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

#include "reflectiveSurfaceMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "photoBioDOM.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    diffuseFraction_(0.0),
    reflectionCoef_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const reflectiveSurfaceMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    diffuseFraction_(ptf.diffuseFraction_),
    reflectionCoef_(ptf.reflectionCoef_)
{}


Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction"))),
        reflectionCoef_(readScalar(dict.lookup("reflectionCoef")))
{


    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
           
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
              
    }
    else
    {
        // No value given. Restart as fixedValue b.c.

        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
        
    }
    
//    Info << "\n n1  \t" << diffuseFraction_ << endl;

    
}


Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const reflectiveSurfaceMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    diffuseFraction_(ptf.diffuseFraction_),
    reflectionCoef_(ptf.reflectionCoef_)
{}


Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const reflectiveSurfaceMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    diffuseFraction_(ptf.diffuseFraction_),
    reflectionCoef_(ptf.reflectionCoef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

      
    scalarField& Iw = *this;
    
    const photoBioModel& photoBio = db().lookupObject<photoBioModel>("photoBioProperties");

    const photoBioDOM& dom(refCast<const photoBioDOM>(photoBio));

    if (dom.nBand() == 0)
    {
        FatalErrorIn
        (
            "Foam::photoBio::"
            "wideBandDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs"
        )   << " a non-grey boundary condition is used with a grey "
            << "absorption model" << nl << exit(FatalError);
    }
    
    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);

    const label patchI = patch().index();
    const  vectorField n = patch().Sf()/patch().magSf();
    
    const  label iBand = dom.IRay(rayId).iBand();
    const  label nAngle = dom.nAngle();     
    const  label nTheta = dom.nTheta();    
    const  label nPhi = dom.nPhi();    
    const  vector& bdRayDir = dom.IRay(rayId).d();
    const  scalar bdOmega  = dom.IRay(rayId).omega();
    
    const scalar deltaPhi   =  pi /(2.0*nPhi);
    const scalar deltaTheta =  pi  /nTheta;
    label  npPhi  = dom.nPixelPhi();
    label  npTheta = dom.nPixelTheta();
    const scalar bdRayPhi  = dom.IRay(rayId).phi();
    const scalar bdRayTheta  = dom.IRay(rayId).theta();

   // const fvPatchScalarField&  dsFace =  dom.diffusionScatter().boundaryField()[patchI];
    
	if (internalField().mesh().nSolutionD() == 2)    //2D (X & Y)
	{	
		npTheta = 1;
    }
    if (internalField().mesh().nSolutionD() == 1)    //2D (X & Y)
	{	
		npTheta = 1; npPhi =1;
    }
    
  	scalar specular = 0.0;
	scalar diffusive = 0.0;

	
    forAll(Iw, faceI)
    {
         specular = 0.0;
	     diffusive = 0.0;	
	     vector surfNorm = -n[faceI];
         scalar cosA = surfNorm& bdRayDir; 
	
	   if(cosA > 0.0 )    // direction out of the wall   
       {
			if( diffuseFraction_ > 0)
			{	 
				for (label jAngle = 0; jAngle < nAngle; jAngle++)
				{         
				label sweepRayID = jAngle + iBand*nAngle;
				vector sweepDir = dom.IRay(sweepRayID).d();
				vector sweepdAve = dom.IRay(sweepRayID).dAve();   
				
				scalar cosB = surfNorm& sweepDir; 
			
				if(  cosB > 0.0 )      // direction out of the wall   
				{ 
				    vector reflecIncidentDir = sweepDir - 2*cosB*surfNorm;		
					label reflecIncidentRay = -1;
                    dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                    const scalarField&  reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];     		
       
					diffusive = diffusive + reflecFace[faceI]*mag(surfNorm & sweepdAve) ;  
				}
				}
			}


		   for (label i = 1; i <= npTheta; i++)
           {
			scalar pxRayTheta = bdRayTheta - 0.5*deltaTheta + 0.5*(2*i -1)*deltaTheta/npTheta;
				    
			for (label j = 1; j <= npPhi; j++)
			{
				scalar pxRayPhi = bdRayPhi - 0.5*deltaPhi + 0.5*(2*j -1)*deltaPhi/npPhi;
					
				scalar sinTheta = Foam::sin(pxRayTheta);
				scalar cosTheta = Foam::cos(pxRayTheta);
				scalar sinPhi = Foam::sin(pxRayPhi);
				scalar cosPhi = Foam::cos(pxRayPhi);
					
				vector  pixelDir = vector(sinTheta*cosPhi , sinTheta* sinPhi , cosTheta);
					
				scalar cosB =  pixelDir & surfNorm;
					
				if ( cosB > 0.0)   
				{
						
					scalar  pixelOmega = 2.0*sinTheta*Foam::sin(deltaTheta/2.0/npTheta)*deltaPhi/npPhi;

					vector reflecIncidentDir = pixelDir - 2*cosB*surfNorm;		
					
					label reflecIncidentRay = -1;
                    dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                    
                    const  scalarField&  reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];                                	

					specular = specular + reflecFace[faceI]*pixelOmega;   
				
				}
			}
		   }

   //         label  i0 = label(Foam::acos(cosA)/deltaPhi);
            refValue()[faceI] = reflectionCoef_*
					(diffuseFraction_*diffusive/pi/2  + (1.0 - diffuseFraction_)*specular/bdOmega); 
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used
        }
    }

    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
    
}


void Foam::photoBio::reflectiveSurfaceMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("reflectionCoef  ")  << reflectionCoef_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction  ")  << diffuseFraction_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{
    makePatchTypeField
    (
        fvPatchScalarField,
        reflectiveSurfaceMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
