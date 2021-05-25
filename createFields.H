Info<< "Reading field T\n" << endl;

volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.0*(B-B.weightedAverage(mesh.V()))
);


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedVector Umoy ("Umoy", U.weightedAverage(mesh.V()));

Info<< "Umoy = " << Umoy.value() << nl << endl;


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar Dmol
(
    transportProperties.lookup("Dmol")
);

#include "createPhi.H"


label  BRefCell = 0;
vector BRefValue(vector::zero);
