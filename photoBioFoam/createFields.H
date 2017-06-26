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

autoPtr<photoBio::photoBioModel> photoBioM
    (
       photoBio::photoBioModel::New(IDefault)
    );
