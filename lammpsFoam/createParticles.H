// Last updated Time-stamp: <2009-03-26 16:00:21 xiao>

Info << "Constructing particles and cloud\n" << endl;

// In MD section, alpha is only kept as a reference. No space
// allocated. Actual field was kept in fluid section.
enhancedCloud cloud
(
    U,
    p,
    Ua,
    Ub,
    nub,
    alpha,
    cloudProperties,
    transportProperties,
    diffusionBandWidth,
    diffusionSteps
);

// Also update/create MD-related fields.
beta = 1.0 - alpha;

volScalarField dragCoef =
    alpha*dimensionedScalar("dum", dimensionSet(1, -3, -1, 0, 0), 0.0);

// Lift (place holder only. Changed in fluid loop)
volVectorField liftCoeff = Cl*beta*rhob*(Ub ^ fvc::curl(U));

Info << "Finished creating particles." << endl;
