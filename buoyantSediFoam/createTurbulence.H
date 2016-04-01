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

    surfaceScalarField betaf("betaf", fvc::interpolate(beta));

    surfaceScalarField betaPhib("betaPhib", betaf*phib);
    autoPtr<PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel> >
    continuousPhaseTurbulence
    (
        PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
        (
            beta,
            Ub,
            betaPhib,
            phib,
            continuousPhaseTransport
        )
    );
