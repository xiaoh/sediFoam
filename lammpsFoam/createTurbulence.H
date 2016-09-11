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

    // Adding the DNS force
    // This aims to simulate the turbulent flow in periodic box
    
    Switch addDNSForce(transportProperties.lookupOrDefault<Switch>("addDNSForce", false));

    if (addDNSForce)
    {
        Info<< "NOTICE!! Adding DNS force in present simulations." << endl;
    }

    volVectorField turbulenceForce
    (
        Ub/dimensionedScalar("dt", dimTime, runTime.deltaTValue())
    );

    turbulenceForce.internalField() *= 0.0;

#ifdef RANDOM_TURB
    Kmesh K(mesh);
    UOprocess forceGen(K, runTime.deltaTValue(), transportProperties);
#endif
