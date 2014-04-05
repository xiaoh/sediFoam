    IOdictionary turbulenceProperties
    (
        IOobject
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Switch turbulence
    (
        turbulenceProperties.lookup("turbulence")
    );

    dictionary kEpsilonCoeffs
    (
        turbulenceProperties.subDict("kEpsilonCoeffs")
    );


    dimensionedScalar Cmu
    (
        kEpsilonCoeffs.lookup("Cmu")
    );

    dimensionedScalar C1
    (
        kEpsilonCoeffs.lookup("C1")
    );

    dimensionedScalar C2
    (
        kEpsilonCoeffs.lookup("C2")
    );

    dimensionedScalar alphak
    (
        kEpsilonCoeffs.lookup("alphak")
    );

    dimensionedScalar alphaEps
    (
        kEpsilonCoeffs.lookup("alphaEps")
    );

    dictionary wallFunctionCoeffs
    (
        turbulenceProperties.subDict("wallFunctionCoeffs")
    );

    dimensionedScalar kappa
    (
        wallFunctionCoeffs.lookup("kappa")
    );

    dimensionedScalar E
    (
        wallFunctionCoeffs.lookup("E")
    );

    nearWallDist y(mesh);


    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field nutb\n" << endl;
    volScalarField nutb
    (
        IOobject
        (
            "nutb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu*sqr(k)/epsilon
    );

    /*  Info<< "Calculating field nuEffa\n" << endl;
    volScalarField nuEffa
    (
        IOobject
        (
            "nuEffa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(Ct)*nutb + nua
        ); */

    Info<< "Calculating field nuEffb\n" << endl;
    volScalarField nuEffb
    (
        IOobject
        (
            "nuEffb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nutb + nub
    );
