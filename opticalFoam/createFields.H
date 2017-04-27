

            volScalarField  IDefault
            (
                IOobject
                (
                    "IDefault",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
               mesh
            );




    autoPtr<optical::opticalModel> opticalM
    (
       optical::opticalModel::New(IDefault)
    );

