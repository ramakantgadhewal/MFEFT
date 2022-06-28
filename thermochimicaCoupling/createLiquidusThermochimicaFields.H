    volScalarField PureCondPhase
    (
        IOobject
        (
            "PureCondPhase",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField FLiNaK_Liquid
    (
        IOobject
        (
            "FLiNaK_Liquid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField K_liquid
    (
        IOobject
        (
            "K_liquid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField Cs_Liquid
    (
        IOobject
        (
            "Cs_Liquid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );
