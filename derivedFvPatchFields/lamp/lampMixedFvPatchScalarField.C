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

    #include "lampMixedFvPatchScalarField.H"
    #include "addToRunTimeSelectionTable.H"
    #include "fvPatchFieldMapper.H"
    #include "volFields.H"

    #include "opticalDOM.H"
    #include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    mixedFvPatchScalarField(p, iF),
    I0_(0.0), 
    nBands_(1),
    lampRadius_(0.0),
    reflectionOnSurface_(false),
    reflectionCoef_(0.0),
    diffuseFraction_(0.0)
    //   opticalBandDist_(null),
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::optical::lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const lampMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    lampRadius_(ptf.lampRadius_),
    reflectionOnSurface_(ptf.reflectionOnSurface_),
    reflectionCoef_(ptf.reflectionCoef_),
    diffuseFraction_(ptf.diffuseFraction_),
    opticalBandDist_(ptf.opticalBandDist_)
{}


Foam::optical::lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    mixedFvPatchScalarField(p, iF),
    I0_(readScalar(dict.lookup("irradiation"))),
    nBands_(readLabel(dict.lookup("nBands"))), 
    lampRadius_(readScalar(dict.lookup("lampRadius"))),
    reflectionCoef_(readScalar(dict.lookup("reflectionCoef"))),
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction")))
{
	
    dict.lookup("reflectionOnSurface") >>reflectionOnSurface_;  
   
    opticalBandDist_.setSize(nBands_);   
    dict.lookup("opticalBandDist") >> opticalBandDist_;
   
	
    if (dict.found("value"))
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
    
}


Foam::optical::lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const lampMixedFvPatchScalarField& ptf
)
    :
    mixedFvPatchScalarField(ptf),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    lampRadius_(ptf.lampRadius_),
    reflectionOnSurface_(ptf.reflectionOnSurface_),
    reflectionCoef_(ptf.reflectionCoef_),
    diffuseFraction_(ptf.diffuseFraction_),
    opticalBandDist_(ptf.opticalBandDist_)
{}


Foam::optical::lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const lampMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    mixedFvPatchScalarField(ptf, iF),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    lampRadius_(ptf.lampRadius_),
    reflectionOnSurface_(ptf.reflectionOnSurface_),
    reflectionCoef_(ptf.reflectionCoef_),
    diffuseFraction_(ptf.diffuseFraction_),
    opticalBandDist_(ptf.opticalBandDist_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::lampMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

     	
    scalarField& Iw = *this;
   
    const opticalModel& optical = db().lookupObject<opticalModel>("opticalProperties");

    const opticalDOM& dom(refCast<const opticalDOM>(optical));
    
    if (dom.nBand() == 0 )  // got rid of the hex mesh requirement
    {
        FatalErrorIn
            (
                "Foam::optical::"
                "wideBandDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs"
            )   << " a non-grey boundary condition is used with a grey "
                << "absorption model" << nl << exit(FatalError);
    }
    const label patchI = patch().index();
    const  vectorField n = patch().Sf()/patch().magSf();
    const  scalarField sfSize = patch().magSf();
    
    label rayId = -1;
    dom.setRayId(dimensionedInternalField().name(), rayId);
    
    const  label nAngle = dom.nAngle();    

    const  label iBand = dom.IRay(rayId).iBand();
    const  vector& bdRayDir = dom.IRay(rayId).d();

    const  scalar bdOmega = dom.IRay(rayId).omega();
    const  scalar bdRayPhi  = dom.IRay(rayId).phi();
    const  scalar bdRayTheta = dom.IRay(rayId).theta();

    const  label nTheta = dom.nTheta();    
    const  label nPhi = dom.nPhi();    
    const  scalar deltaPhi   =  pi /(2.0*nPhi);
    const  scalar deltaTheta =  pi  /nTheta;

    label  npPhi  = dom.NumPixelPhi();    
    label  npTheta = dom.NumPixelTheta();  

    const scalar nOwn_ = 1.00;        // refractive index of the internal fluid solution (air as default) 
    const scalar nNbg_ = 1.5168;      // refractive index of the lamp material
    const scalar nRatio = nOwn_/  nNbg_;
    const scalar height = 0.001;     
     
    angleDist.setSize(nPhi); 
    scalar sumP =0;      
    for(label i = 0;i < nPhi;i++)
    {	
        scalar alpha = (2*i+1)/2*deltaPhi ;
        angleDist[i] = alpha + Foam::sin(alpha)*Foam::cos(alpha)- sumP ;  
        sumP = sumP + angleDist[i];
        angleDist[i]  = angleDist[i] /pi;
    }
    angleDist[0] = 2*angleDist[0];
   
    scalar specularReflection = 0.0;
    scalar reflectionRefraction = 0.0;
    scalar  diffuseRefractionInside = 0.0;	
    autoPtr<scalarField>   diffusiveReflection;	   
    diffusiveReflection.set(new scalarField(patch().size(),0.0)); 
     
    if(reflectionOnSurface_)
    {     
	
        forAll(Iw, faceI)
        {
            vector surfNorm = -n[faceI] ;
       
            for (label jAngle = 0; jAngle < nAngle; jAngle++)
            {         
                label sweepRayID = jAngle + (iBand)*nAngle;
                vector sweepDir = dom.IRay(sweepRayID).d();
                vector sweepdAve = dom.IRay(sweepRayID).dAve();   
                scalar sweepdOmega = dom.IRay(sweepRayID).omega();   
     
                scalar cosA = surfNorm & sweepDir; 
     
                if(  cosA < 0.0  && cosA >= -1)      // direction into  the wall   
                {
                    const scalarField&  reflecFace = dom.IRay(sweepRayID).I().boundaryField()[patchI];            
                 
                    scalar sinA = Foam::sqrt(1-cosA*cosA);
                    scalar sinB = sinA*nRatio ; 
                    scalar cosB =-Foam::sqrt(1-sinB*sinB); 
    
                    scalar r1 = (nNbg_*cosB - nOwn_*cosA ) / (nNbg_*cosB + nOwn_*cosA );
                    scalar r2 = (nNbg_*cosA - nOwn_*cosB ) / (nNbg_*cosA + nOwn_*cosB );
     
                    scalar R = 0.5*(r1*r1 + r2*r2);
     
                    diffusiveReflection()[faceI] = diffusiveReflection()[faceI]+ reflecFace[faceI]* R *mag( n[faceI] & sweepdAve); 
                    diffuseRefractionInside = diffuseRefractionInside + reflecFace[faceI] * ( 1.0 - R )* sfSize[faceI] * sweepdOmega ;  //change from
                }
            }
        }
    }
    
     
     
    forAll(Iw, faceI)
    {
        vector surfNorm = -n[faceI] ;
     
        scalar cosA = surfNorm & bdRayDir; 
        label  i0 = label(Foam::acos(cosA)/deltaPhi);
		
        if(cosA > 0)
        {  
            reflectionRefraction = 0.0;
		 
            if(reflectionOnSurface_)
            {
                specularReflection = 0.0;

                scalar pxRayTheta = bdRayTheta;
     
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
     
                        const scalarField&  reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];                                	
     
     
                        if( cosB*cosB > 1-1/(nRatio*nRatio) )      // whether  total reflection or not  
                        { 
     
                            vector refracIncidentDir = (pixelDir - cosB * surfNorm )*nRatio 
                                + Foam::sqrt(1 -(nRatio*nRatio)*(1- cosB*cosB))*surfNorm ;
     
                            scalar cosA = surfNorm& refracIncidentDir; // / mag(d);                //    /mag(n[faceI])
     
                            scalar r1 = (nNbg_*cosB - nOwn_*cosA )/ (nNbg_*cosB + nOwn_*cosA );
                            scalar r2 = (nNbg_*cosA - nOwn_*cosB )/ (nNbg_*cosA + nOwn_*cosB );
     
                            scalar R =  0.5*(r1*r1 + r2*r2);
     
                            specularReflection = specularReflection + reflecFace[faceI]*R*pixelOmega;
     
                        }   
                        else
                        {
                            specularReflection = specularReflection + reflecFace[faceI]*pixelOmega;   
                        }
                    }
                }
     
                reflectionRefraction =  diffuseRefractionInside/(2*pi*lampRadius_*height)/pi/2            
                    + (diffuseFraction_*diffusiveReflection()[faceI] /pi/2 
                    + (1.0 - diffuseFraction_)*specularReflection/bdOmega) ;
     
            }  
          
     
            refValue()[faceI] = reflectionCoef_*reflectionRefraction  + I0_*opticalBandDist_[iBand]*angleDist[i0]/bdOmega;  ///pi/2;  //
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

void Foam::optical::lampMixedFvPatchScalarField::dirToAngle
(
    const vector& dir,
    scalar&	tPhi,
    scalar& tTheta
) const
{

    tTheta = Foam::acos(dir.z()/mag(dir));


    if(dir.x() != 0 )
    { 
        tPhi = Foam::atan(dir.y()/dir.x());
        if(dir.x() < 0 && dir.y() >= 0 ) tPhi = tPhi + pi;
        if(dir.x() < 0 && dir.y() < 0 ) tPhi = tPhi + pi;
        if(dir.x() > 0 && dir.y() < 0 ) tPhi = tPhi + 2*pi;
    }
    else
    {
        if(dir.y() > 0 ) tPhi = pi/2.0;
        if(dir.y() < 0 ) tPhi = 3.0*pi/2.0;;
    }
	
}

void Foam::optical::lampMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("irradiation") << I0_ << token::END_STATEMENT << nl;
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("lampRadius") << lampRadius_ << token::END_STATEMENT << nl;
    os.writeKeyword("reflectionOnSurface") << reflectionOnSurface_ << token::END_STATEMENT << nl;
    os.writeKeyword("reflectionCoef") << reflectionCoef_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("opticalBandDist") << opticalBandDist_ << token::END_STATEMENT << nl;                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        makePatchTypeField
        (
            fvPatchScalarField,
            lampMixedFvPatchScalarField
        );
    }
}


// ************************************************************************* //
