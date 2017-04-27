    // Initialise field pointer lists
    
    PtrList<volScalarField> I_cal(calRegions.size());
    
    PtrList<optical::opticalModel> opticalModels(calRegions.size());

    // Populate solid field pointer lists
    forAll(calRegions, i)
    {
        Info<< "*** Reading  properties for region "
            << calRegions[i].name() << nl << endl;

        I_cal.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "IDefault",
                    runTime.timeName(),
                    calRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                calRegions[i]
            )
        );

        opticalModels.set
        (
            i,
            optical::opticalModel::New(I_cal[i])
        );
    }
    
    
    
