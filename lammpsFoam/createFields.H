    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1) - alpha
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Ub\n" << endl;
    volVectorField Ub
    (
        IOobject
        (
            "Ub",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Setting up place-holder for Ua\n" << endl;
    volVectorField Ua
    (
        IOobject
        (
            "Ua",
            runTime.timeName(),
            mesh,
            // IOobject::READ_IF_PRESENT,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        Ub - Ub
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha*Ua + beta*Ub
    );

    // Ua*alpha (sediment flux)
    volVectorField Uc
    (
        IOobject
        (
            "Uc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha*Ua
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary cloudProperties
    (
        IOobject
        (
            "cloudProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // default diffusion band width if not defined in dictionary
    scalar diffusionBandWidth = 0.006;
    if (cloudProperties.found("diffusionBandWidth"))
    {
        diffusionBandWidth =
            readScalar
            (
                cloudProperties.lookup("diffusionBandWidth")
            );
    }
    else
    {
        Info<< "*** Diffusion band width set to 0.006 by default"
            << endl;
    }

    // default max alpha if not defined in dictionary
    scalar maxPossibleAlpha = 0.70;
    if (cloudProperties.found("maxPossibleAlpha"))
    {
        maxPossibleAlpha =
            readScalar
            (
                cloudProperties.lookup("maxPossibleAlpha")
            );
    }
    else
    {
        Info<< "*** Max possible alpha set to 0.70 by default"
            << endl;
    }

    // default diffusion steps if not defined in dictionary
    scalar diffusionSteps = 6;
    if (cloudProperties.found("diffusionSteps"))
    {
        diffusionSteps =
            readScalar
            (
                cloudProperties.lookup("diffusionSteps")
            );
    }
    else
    {
        Info<< "*** Diffusion steps set to 6 by default"
            << endl;
    }

    dimensionedScalar rhob
    (
        transportProperties.lookup("rhob")
    );

    dimensionedScalar nub
    (
        transportProperties.lookup("nub")
    );

    dimensionedScalar Cvm
    (
        transportProperties.lookup("Cvm")
    );

    dimensionedScalar Cl
    (
        transportProperties.lookup("Cl")
    );

    dimensionedScalar Ct
    (
        transportProperties.lookup("Ct")
    );

#   include "createPhia.H"
#   include "createPhib.H"

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh
        ),
        fvc::interpolate(alpha)*phia
      + fvc::interpolate(beta)*phib
    );


#   include "createTurbulence.H"

    Info<< "Calculating field DDtUa and DDtUb\n" << endl;

    volVectorField DDtUa =
        fvc::ddt(Ua)
      + fvc::div(phia, Ua)
      - fvc::div(phia)*Ua;

    volVectorField DDtUb =
        fvc::ddt(Ub)
      + fvc::div(phib, Ub)
      - fvc::div(phib)*Ub;


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());

    // Adding pressure gradient
    chPressureGrad::initPressureGrad(transportProperties);

    chPressureGrad gradP(Ub,alpha,"gradP","LAMMPS");

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
