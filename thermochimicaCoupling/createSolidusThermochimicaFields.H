    volScalarField FLiNaK_Solid
    (
    IOobject
    (
        "FLiNaK_Solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("", dimless, 0)
    );

    volScalarField solid1
    (
        IOobject
        (
            "Phase_LiF_NaF_KF_CsF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField NaF_KF
    (
        IOobject
        (
            "Phase_NaF_KF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField LiF_NaF
    (
        IOobject
        (
            "Phase_LiF_NaF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField LiF
    (
        IOobject
        (
            "LiF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField NaF
    (
        IOobject
        (
            "NaF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField KF
    (
        IOobject
        (
            "KF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField CsF
    (
        IOobject
        (
            "CsF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField LaK3F6
    (
        IOobject
        (
            "LaK3F6_solid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField CeF3
    (
        IOobject
        (
            "CeF3_solid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField CsF2Li_solids
    (
        IOobject
        (
            "CsF2Li_solids",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField FLi_solids
    (
        IOobject
        (
            "FLi_solids",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField CsF_Solids
    (
        IOobject
        (
            "CsF_Solids",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );
